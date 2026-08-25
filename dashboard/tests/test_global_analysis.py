from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from rare_oncology_dashboard.catalog import (
    ANALYSIS_KIND,
    THERAPY_KIND,
    ExecutionSnapshot,
)
from rare_oncology_dashboard.global_analysis import (
    build_global_overview,
    convergence_matrix,
    load_family_catalog,
)


def snapshot(
    execution_id: int,
    kind: str,
    gene: str,
    created_at: str,
    result: dict,
    *,
    parameters: dict | None = None,
) -> ExecutionSnapshot:
    return ExecutionSnapshot(
        database_path=Path(f"/{gene.lower()}.sqlite3"),
        execution_id=execution_id,
        kind=kind,
        gene=gene,
        parameter_hash=f"hash-{execution_id}",
        parameters=parameters or {},
        result=result,
        created_at=created_at,
        last_accessed_at=created_at,
        cache_hit_count=0,
        schema_version=1,
    )


def analysis_result(gene: str, change: str, level: str) -> dict:
    return {
        "gene": gene,
        "rare_variant_count": 4,
        "cancer_observation_count": 9,
        "shared_mutations": [
            {
                "protein_change": change,
                "diseases": [f"{gene} rare disease"],
                "orphanet_codes": ["ORPHA:1"],
                "match_level": level,
                "cancer_patient_count": 7,
                "tcga_projects": ["TCGA-LUAD"],
            }
        ],
        "source_statuses": [
            {"source": "ClinVar", "status": "ok"},
            {"source": "GDC", "status": "partial"},
        ],
    }


def signature(result: dict) -> dict:
    record = result["shared_mutations"][0]
    return {
        "gene": result["gene"],
        "shared_mutations": [
            {
                "protein_change": record["protein_change"],
                "diseases": record["diseases"],
                "orphanet_codes": record["orphanet_codes"],
                "match_level": record["match_level"],
            }
        ],
    }


class GlobalAnalysisTests(unittest.TestCase):
    def test_uses_newest_analysis_and_only_exactly_matched_therapy(self) -> None:
        old_egfr = analysis_result("EGFR", "T790M", "candidate_protein_match")
        new_egfr = analysis_result(
            "EGFR", "L858R", "confirmed_genomic_and_protein"
        )
        new_egfr["shared_mutations"][0]["diseases"] = [
            "Hereditary cancer-predisposing syndrome",
            "EGFR rare disease",
        ]
        acvr1 = analysis_result("ACVR1", "R206H", "candidate_protein_match")
        analyses = [
            snapshot(3, ANALYSIS_KIND, "EGFR", "2026-03-01T00:00:00+00:00", new_egfr),
            snapshot(2, ANALYSIS_KIND, "ACVR1", "2026-02-01T00:00:00+00:00", acvr1),
            snapshot(1, ANALYSIS_KIND, "EGFR", "2026-01-01T00:00:00+00:00", old_egfr),
        ]
        incompatible = snapshot(
            12,
            THERAPY_KIND,
            "EGFR",
            "2026-04-01T00:00:00+00:00",
            {
                "variant_therapy_evidence": [{"source": "CIViC"}],
                "regulatory_therapy": [],
                "trial_and_target_discovery": [],
            },
            parameters={"arguments": {"analysis_signature": signature(old_egfr)}},
        )
        compatible = snapshot(
            11,
            THERAPY_KIND,
            "EGFR",
            "2026-03-02T00:00:00+00:00",
            {
                "variant_therapy_evidence": [{"source": "CIViC"}],
                "regulatory_therapy": [
                    {"source": "FDA", "therapy_status": "regulatory_approved"}
                ],
                "trial_and_target_discovery": [
                    {"source": "CTG", "therapy_status": "clinical_trial"},
                    {"source": "OT", "therapy_status": "mechanism_target"},
                ],
            },
            parameters={"arguments": {"analysis_signature": signature(new_egfr)}},
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "receptor_tyrosine_kinases.txt").write_text("EGFR\n")
            (root / "tgfb_receptors.txt").write_text("ACVR1\n")
            (root / "kinase_linked_cytokine_receptors.txt").write_text("")
            families = load_family_catalog(root)

        overview = build_global_overview(
            analyses, [incompatible, compatible], families
        )

        self.assertEqual(overview.active_gene_count, 2)
        self.assertEqual(overview.analysis_snapshot_count, 3)
        self.assertEqual(overview.shared_variant_count, 2)
        self.assertEqual(overview.confirmed_shared_count, 1)
        self.assertEqual(overview.matched_therapy_count, 1)
        egfr = next(row for row in overview.gene_summaries if row.gene == "EGFR")
        self.assertEqual(egfr.therapy_execution_id, 11)
        self.assertEqual(egfr.variant_therapy_count, 1)
        self.assertEqual(egfr.actionability_tier, "Variant-specific evidence")
        self.assertEqual(egfr.hereditary_cancer_disease_count, 1)
        self.assertEqual(egfr.other_inherited_disease_count, 1)
        acvr1_row = next(row for row in overview.gene_summaries if row.gene == "ACVR1")
        self.assertEqual(acvr1_row.protein_family, "TGF-beta type I receptors")
        self.assertEqual(acvr1_row.actionability_tier, "No therapy evidence")
        self.assertIn(
            ("EGFR rare disease", "EGFR/ERBB family", 1),
            overview.disease_family_counts,
        )
        self.assertIn(
            (
                "Hereditary cancer-predisposing syndrome",
                "EGFR/ERBB family",
                1,
            ),
            overview.hereditary_cancer_family_counts,
        )
        self.assertIn(
            ("EGFR rare disease", "EGFR/ERBB family", 1),
            overview.other_inherited_disease_family_counts,
        )

    def test_convergence_matrix_is_ranked_and_bounded(self) -> None:
        matrix = convergence_matrix(
            [
                ("Disease B", "Family 1", 2),
                ("Disease A", "Family 1", 4),
                ("Disease B", "Family 2", 3),
                ("Disease C", "Family 2", 1),
            ],
            top_rows=2,
        )

        self.assertEqual(matrix["rows"], ["Disease B", "Disease A"])
        self.assertEqual(matrix["families"], ["Family 1", "Family 2"])
        self.assertEqual(matrix["values"], [[2, 3], [4, 0]])

    def test_rejects_duplicate_curated_family_membership(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "receptor_tyrosine_kinases.txt").write_text("EGFR\n")
            (root / "tgfb_receptors.txt").write_text("EGFR\n")
            (root / "kinase_linked_cytokine_receptors.txt").write_text("")
            with self.assertRaisesRegex(ValueError, "more than one family"):
                load_family_catalog(root)


if __name__ == "__main__":
    unittest.main()
