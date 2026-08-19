from __future__ import annotations

import csv
import hashlib
import json
import sqlite3
import unittest
from pathlib import Path

from rare_oncology.storage import analysis_result_from_dict, therapy_result_from_dict


ROOT = Path(__file__).resolve().parents[1]
BOOTSTRAP = ROOT / "bootstrap"
DATABASE = BOOTSTRAP / "rare_oncology_demo.sqlite3"
GENES = {"EGFR", "TGFBR1", "ACVR1", "PDGFRB"}


class BootstrapBundleTests(unittest.TestCase):
    def test_manifest_files_are_small_and_match_sha256(self) -> None:
        manifest = json.loads((BOOTSTRAP / "manifest.json").read_text())
        self.assertEqual(set(manifest["genes"]), GENES)
        descriptors = [manifest["database"], *manifest["regulatory_tables"]]
        for descriptor in descriptors:
            path = BOOTSTRAP / descriptor["path"]
            self.assertTrue(path.is_file())
            self.assertEqual(path.stat().st_size, descriptor["size_bytes"])
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(),
                descriptor["sha256"],
            )
        self.assertLess(DATABASE.stat().st_size, 10 * 1024 * 1024)

    def test_seed_contains_rehydratable_analysis_and_compatible_therapy(self) -> None:
        with sqlite3.connect(f"file:{DATABASE}?mode=ro", uri=True) as connection:
            connection.row_factory = sqlite3.Row
            rows = connection.execute(
                """
                SELECT kind, gene, parameters_json, result_json
                FROM execution_cache ORDER BY gene, kind
                """
            ).fetchall()
        self.assertEqual(len(rows), 8)
        self.assertEqual({str(row["gene"]) for row in rows}, GENES)
        by_gene = {
            gene: {str(row["kind"]): row for row in rows if row["gene"] == gene}
            for gene in GENES
        }
        for gene, records in by_gene.items():
            analysis_payload = json.loads(records["mutation_analysis"]["result_json"])
            therapy_payload = json.loads(records["therapy_analysis"]["result_json"])
            analysis = analysis_result_from_dict(analysis_payload)
            therapy = therapy_result_from_dict(therapy_payload)
            self.assertEqual(analysis.gene, gene)
            self.assertEqual(therapy.gene, gene)
            self.assertIn("Bootstrap snapshot", analysis.caveat)
            self.assertTrue(analysis.shared_mutations)
            self.assertTrue(analysis.mutation_plot_data)
            self.assertLessEqual(
                {record.source.casefold() for record in analysis.evidence_records},
                {"alphamissense", "gnomad", "panelapp"},
            )
            parameters = json.loads(records["therapy_analysis"]["parameters_json"])
            self.assertEqual(
                parameters["arguments"]["analysis_signature"],
                _analysis_signature(analysis_payload),
            )

    def test_seed_contains_available_string_networks_and_no_builder_path(self) -> None:
        payload = DATABASE.read_bytes()
        self.assertFalse(b"/Users/" in payload, "Seed contains a machine-local path")
        with sqlite3.connect(f"file:{DATABASE}?mode=ro", uri=True) as connection:
            rows = connection.execute(
                "SELECT gene, result_json FROM string_network_cache ORDER BY gene"
            ).fetchall()
        self.assertEqual({row[0] for row in rows}, {"ACVR1", "PDGFRB", "TGFBR1"})
        for gene, result_json in rows:
            result = json.loads(result_json)
            self.assertEqual(result["gene"], gene)
            self.assertTrue(result["nodes"])
            self.assertTrue(result["edges"])

    def test_normalized_regulatory_tables_are_portable(self) -> None:
        expected = {
            "fda_biomarkers.csv": {"EGFR", "PDGFRB"},
            "ema_labels.csv": {"EGFR"},
        }
        for filename, required_genes in expected.items():
            with (BOOTSTRAP / "regulatory" / filename).open(
                encoding="utf-8-sig", newline=""
            ) as handle:
                reader = csv.DictReader(handle)
                self.assertIn("source_url", reader.fieldnames or [])
                self.assertIn("retrieved_at", reader.fieldnames or [])
                genes = {str(row["gene"]) for row in reader}
            self.assertLessEqual(required_genes, genes)


def _analysis_signature(result: dict) -> dict:
    return {
        "gene": str(result.get("gene") or "").upper(),
        "shared_mutations": [
            {
                "protein_change": record.get("protein_change"),
                "diseases": list(record.get("diseases") or []),
                "orphanet_codes": list(record.get("orphanet_codes") or []),
                "match_level": record.get("match_level"),
            }
            for record in result.get("shared_mutations", [])
        ],
    }


if __name__ == "__main__":
    unittest.main()
