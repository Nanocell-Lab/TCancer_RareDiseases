from __future__ import annotations

import json
import tempfile
import unittest
from contextlib import redirect_stdout
from datetime import datetime, timedelta, timezone
from io import StringIO
from pathlib import Path

from rare_oncology.global_analysis import (
    build_global_family_analysis,
    load_family_catalog,
    write_global_report,
)
from rare_oncology.global_cli import main
from rare_oncology.global_plots import generate_global_family_plots
from rare_oncology.models import (
    AnalysisResult,
    SharedMutation,
    SourceStatus,
    TherapyAssociation,
    TherapyResult,
)
from rare_oncology.storage import ANALYSIS_KIND, THERAPY_KIND, ExecutionDatabase


ROOT = Path(__file__).resolve().parents[1]


class Clock:
    def __init__(self) -> None:
        self.value = datetime(2026, 8, 17, tzinfo=timezone.utc)

    def __call__(self) -> datetime:
        current = self.value
        self.value += timedelta(seconds=1)
        return current


def analysis_result(
    gene: str = "EGFR", change: str = "L858R", *, patients: int = 4
) -> AnalysisResult:
    shared = SharedMutation(
        gene=gene,
        protein_change=change,
        diseases=(f"Rare {gene} disorder",),
        orphanet_codes=("999",),
        clinvar_ids=("101",),
        clinvar_accessions=("VCV101",),
        germline_classifications=("Pathogenic",),
        canonical_spdis=("NC_000007.14:1:A:T",),
        tcga_projects=("TCGA-LUAD",),
        cancer_patient_count=patients,
        cancer_event_count=patients,
        genomic_match_count=1,
        genomic_confirmations=("GRCh38:7:2:A>T",),
        ssm_ids=("ssm-1",),
        match_level="confirmed_genomic_and_protein",
    )
    return AnalysisResult(
        gene=gene,
        generated_at="2026-08-17T00:00:00+00:00",
        gdc_program="TCGA",
        gdc_projects=("TCGA-LUAD",),
        min_cancer_patients=1,
        excluded_conditions=(),
        clinvar_record_limit=None,
        allowed_germline_classifications=("Pathogenic",),
        gdc_data_release="46",
        rare_variant_count=3,
        cancer_observation_count=patients,
        shared_mutation_count=1,
        shared_mutations=(shared,),
        mutation_plot_data=(),
        source_statuses=(
            SourceStatus("ClinVar", "ok", "variant_assertion", 3),
            SourceStatus("GDC", "ok", "somatic_observation", patients),
        ),
    )


def analysis_signature(result: AnalysisResult) -> dict[str, object]:
    return {
        "gene": result.gene,
        "shared_mutations": [
            {
                "protein_change": item.protein_change,
                "diseases": list(item.diseases),
                "orphanet_codes": list(item.orphanet_codes),
                "match_level": item.match_level,
            }
            for item in result.shared_mutations
        ],
    }


def therapy_result(gene: str = "EGFR", change: str = "L858R") -> TherapyResult:
    common = {
        "gene": gene,
        "therapy_components": (),
        "retrieved_at": "2026-08-17T00:00:00+00:00",
    }
    return TherapyResult(
        gene=gene,
        generated_at="2026-08-17T00:00:00+00:00",
        protein_changes=(change,),
        variant_therapy_evidence=(
            TherapyAssociation(
                layer="variant_therapy_evidence",
                source="CIViC",
                source_record_id="civic-1",
                association_scope="candidate_protein_variant",
                therapy_name="osimertinib",
                therapy_status="accepted_predictive_evidence",
                protein_change=change,
                **common,
            ),
        ),
        regulatory_therapy=(
            TherapyAssociation(
                layer="regulatory_therapy",
                source="FDA",
                source_record_id="fda-1",
                association_scope="label_biomarker",
                therapy_name="osimertinib",
                therapy_status="regulatory_approved",
                **common,
            ),
        ),
        trial_and_target_discovery=(
            TherapyAssociation(
                layer="trial_and_target_discovery",
                source="ClinicalTrials.gov",
                source_record_id="NCT1",
                association_scope="trial_variant_mention",
                therapy_name="study intervention",
                therapy_status="clinical_trial",
                **common,
            ),
            TherapyAssociation(
                layer="trial_and_target_discovery",
                source="Open Targets",
                source_record_id="CHEMBL1",
                association_scope="gene_target",
                therapy_name="target compound",
                therapy_status="mechanism_target",
                **common,
            ),
        ),
        source_statuses=(),
    )


class GlobalFamilyTests(unittest.TestCase):
    def test_family_catalog_classifies_curated_receptors(self):
        catalog = load_family_catalog(ROOT / "gene_lists")
        self.assertEqual(len(catalog), 114)
        self.assertEqual(catalog["EGFR"].protein_family, "EGFR/ERBB family")
        self.assertEqual(catalog["TGFBR1"].protein_family, "TGF-beta type I receptors")
        self.assertEqual(catalog["IL6ST"].protein_family, "IL-6 receptor family")

    def test_incremental_index_matches_therapy_to_exact_analysis(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "runs.sqlite3"
            index = root / "global.sqlite3"
            molecular = analysis_result()
            with ExecutionDatabase(source, clock=Clock()) as database:
                database.save(ANALYSIS_KIND, "EGFR", {"run": 1}, molecular)
                database.save(
                    THERAPY_KIND,
                    "EGFR",
                    {"arguments": {"analysis_signature": analysis_signature(molecular)}},
                    therapy_result(),
                )
            first = build_global_family_analysis(
                (source,), index_database=index, gene_list_dir=ROOT / "gene_lists"
            )
            second = build_global_family_analysis(
                (source,), index_database=index, gene_list_dir=ROOT / "gene_lists"
            )
        self.assertEqual(first.indexed_new_snapshots, 2)
        self.assertEqual(second.indexed_new_snapshots, 0)
        self.assertEqual(first.active_analysis_count, 1)
        self.assertEqual(first.active_therapy_count, 1)
        self.assertEqual(first.index_error_count, 0)
        self.assertEqual(first.source_error_count, 0)
        self.assertEqual(
            first.family_group_summaries[0].family_group,
            "Receptor tyrosine kinases",
        )
        gene = first.gene_summaries[0]
        self.assertEqual(gene.protein_family, "EGFR/ERBB family")
        self.assertEqual(gene.confirmed_shared_count, 1)
        self.assertEqual(gene.variant_specific_therapy_count, 1)
        self.assertEqual(gene.active_regulatory_therapy_count, 1)
        self.assertEqual(gene.clinical_trial_count, 1)
        self.assertEqual(gene.mechanism_target_count, 1)
        self.assertEqual(gene.actionability_tier, "Active regulatory biomarker")

    def test_new_analysis_does_not_reuse_therapy_from_an_old_signature(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "runs.sqlite3"
            index = root / "global.sqlite3"
            first_analysis = analysis_result()
            clock = Clock()
            with ExecutionDatabase(source, clock=clock) as database:
                database.save(ANALYSIS_KIND, "EGFR", {"run": 1}, first_analysis)
                database.save(
                    THERAPY_KIND,
                    "EGFR",
                    {"arguments": {"analysis_signature": analysis_signature(first_analysis)}},
                    therapy_result(),
                )
            build_global_family_analysis(
                (source,), index_database=index, gene_list_dir=ROOT / "gene_lists"
            )
            newer = analysis_result(change="T790M", patients=8)
            with ExecutionDatabase(source, clock=clock) as database:
                database.save(ANALYSIS_KIND, "EGFR", {"run": 2}, newer)
            report = build_global_family_analysis(
                (source,), index_database=index, gene_list_dir=ROOT / "gene_lists"
            )
        self.assertEqual(report.indexed_new_snapshots, 1)
        self.assertEqual(report.active_therapy_count, 0)
        self.assertEqual(report.unmatched_therapy_gene_count, 1)
        self.assertIsNone(report.gene_summaries[0].therapy_snapshot_key)
        self.assertEqual(report.gene_summaries[0].max_variant_patient_count, 8)

    def test_identical_database_copy_does_not_duplicate_snapshots(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "runs.sqlite3"
            copied = root / "copied.sqlite3"
            with ExecutionDatabase(source) as database:
                database.save(ANALYSIS_KIND, "TGFBR1", {}, analysis_result("TGFBR1", "R487Q"))
            copied.write_bytes(source.read_bytes())
            report = build_global_family_analysis(
                (source, copied),
                index_database=root / "global.sqlite3",
                gene_list_dir=ROOT / "gene_lists",
            )
        self.assertEqual(report.indexed_new_snapshots, 1)
        self.assertEqual(report.active_analysis_count, 1)

    def test_exports_and_figures_are_created(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "runs.sqlite3"
            with ExecutionDatabase(source) as database:
                database.save(ANALYSIS_KIND, "EGFR", {}, analysis_result())
            report = build_global_family_analysis(
                (source,),
                index_database=root / "global.sqlite3",
                gene_list_dir=ROOT / "gene_lists",
            )
            tables = write_global_report(report, root / "exports")
            plots = generate_global_family_plots(
                report, root / "exports" / "plots", formats=("png",)
            )
            payload = json.loads(tables[0].read_text(encoding="utf-8"))
            self.assertEqual(payload["active_analysis_count"], 1)
            self.assertEqual(len(tables), 4)
            self.assertEqual(len(plots), 4)
            self.assertTrue(all(path.stat().st_size > 100 for path in (*tables, *plots)))

    def test_cli_discovers_sources_and_reuses_incremental_index(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source" / "runs.sqlite3"
            source.parent.mkdir()
            with ExecutionDatabase(source) as database:
                database.save(ANALYSIS_KIND, "EGFR", {}, analysis_result())
            arguments = [
                "--results-root", str(root / "source"),
                "--index-database", str(root / "global.sqlite3"),
                "--gene-list-dir", str(ROOT / "gene_lists"),
                "--output-dir", str(root / "exports"),
            ]
            with redirect_stdout(StringIO()) as output:
                self.assertEqual(main(arguments), 0)
                self.assertEqual(main(arguments), 0)
            text = output.getvalue()
        self.assertIn("1 new snapshots indexed", text)
        self.assertIn("0 new snapshots indexed", text)


if __name__ == "__main__":
    unittest.main()
