from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from rare_oncology.regulatory_tables import (
    main,
    parse_ema_json,
    parse_fda_html,
    write_table,
)
from rare_oncology.therapy import RegulatoryTableTherapyProvider, TherapyContext


FDA_FIXTURE = """
<html><body>
<h2>Pharmacogenomic Biomarkers in Drug Labeling</h2>
<table summary="FDA-approved drugs with pharmacogenomic information in their labeling">
  <thead><tr><th>Drug</th><th>Therapeutic Area*</th><th>Biomarker</th><th>Labeling Sections</th></tr></thead>
  <tbody>
    <tr><td><a href="http://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=overview.process&amp;varApplNo=123456"><strong>Exampletinib (1)</strong></a></td><td>Oncology</td><td>EGFR<br>(L858R)</td><td>Indications and Usage, Clinical Studies</td></tr>
    <tr><td><a href="/drug/example"><strong>Unrelated</strong></a></td><td>Metabolism</td><td>CYP2D6</td><td>Clinical Pharmacology</td></tr>
  </tbody>
</table>
<h2>Content current as of:</h2><p><time datetime="2026-03-03T00:00:00Z">03/03/2026</time></p>
</body></html>
"""


def ema_fixture() -> dict[str, object]:
    return {
        "meta": {"total_records": 4, "timestamp": "2026-08-11T18:10:25Z"},
        "data": [
            {
                "category": "Human",
                "name_of_medicine": "Metinib",
                "ema_product_number": "EMEA/H/C/000001",
                "medicine_status": "Authorised",
                "therapeutic_indication": (
                    "Treatment of NSCLC harbouring alterations leading to "
                    "mesenchymal-epithelial transition factor gene exon 14 "
                    "(METex14) skipping."
                ),
                "therapeutic_area_mesh": "Carcinoma, Non-Small-Cell Lung",
                "active_substance": "metinib",
                "conditional_approval": "Yes",
                "orphan_medicine": "No",
                "medicine_url": "https://www.ema.europa.eu/example/metinib",
            },
            {
                "category": "Human",
                "name_of_medicine": "Old-EGFR",
                "ema_product_number": "EMEA/H/C/000002",
                "medicine_status": "Withdrawn",
                "therapeutic_indication": "Treatment of EGFR-positive cancer.",
                "therapeutic_area_mesh": "Neoplasms",
                "active_substance": "old-egfr",
                "medicine_url": "https://www.ema.europa.eu/example/old-egfr",
            },
            {
                "category": "Human",
                "name_of_medicine": "Metabolic drug",
                "ema_product_number": "EMEA/H/C/000003",
                "medicine_status": "Authorised",
                "therapeutic_indication": "Treatment of metabolic disease.",
                "medicine_url": "https://www.ema.europa.eu/example/metabolic",
            },
            {
                "category": "Human",
                "name_of_medicine": "Prior-therapy-only",
                "ema_product_number": "EMEA/H/C/000005",
                "medicine_status": "Authorised",
                "therapeutic_indication": (
                    "Treatment after chemotherapy and anti-EGFR agents."
                ),
                "medicine_url": "https://www.ema.europa.eu/example/prior",
            },
            {
                "category": "Veterinary",
                "name_of_medicine": "Veterinary EGFR",
                "ema_product_number": "EMEA/V/C/000004",
                "medicine_status": "Authorised",
                "therapeutic_indication": "EGFR-positive disease.",
                "medicine_url": "https://www.ema.europa.eu/example/vet",
            },
        ],
    }


class FdaParserTests(unittest.TestCase):
    def test_builds_conservative_label_biomarker_rows(self):
        rows = parse_fda_html(
            FDA_FIXTURE, ["EGFR", "MET"], retrieved_at="2026-08-11T00:00:00Z"
        )
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["gene"], "EGFR")
        self.assertEqual(row["therapy_name"], "Exampletinib")
        self.assertEqual(row["therapy_status"], "label_biomarker")
        self.assertEqual(row["therapy_id"], "123456")
        self.assertEqual(row["source_version"], "2026-03-03T00:00:00Z")
        self.assertTrue(row["source_url"].startswith("https://"))

    def test_missing_official_table_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "pharmacogenomic table"):
            parse_fda_html("<html></html>", ["EGFR"])


class EmaParserTests(unittest.TestCase):
    def test_exact_and_controlled_alias_matches_keep_regulatory_status(self):
        rows = parse_ema_json(
            ema_fixture(),
            ["MET", "EGFR"],
            retrieved_at="2026-08-11T00:00:00Z",
        )
        self.assertEqual([(row["gene"], row["therapy_name"]) for row in rows], [
            ("EGFR", "Old-EGFR"),
            ("MET", "Metinib"),
        ])
        by_gene = {row["gene"]: row for row in rows}
        self.assertEqual(by_gene["MET"]["therapy_status"], "conditional_approval")
        self.assertEqual(by_gene["MET"]["effect"], "biomarker_selected")
        self.assertIn("conditional approval", by_gene["MET"]["evidence_level"])
        self.assertEqual(by_gene["EGFR"]["therapy_status"], "withdrawn")
        self.assertEqual(by_gene["EGFR"]["effect"], "biomarker_selected")
        self.assertEqual(by_gene["EGFR"]["regulatory_source_status"], "Withdrawn")

    def test_substring_does_not_turn_metabolic_into_met(self):
        rows = parse_ema_json(ema_fixture(), ["MET"])
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["therapy_name"], "Metinib")

    def test_prior_anti_gene_therapy_is_not_a_biomarker_association(self):
        rows = parse_ema_json(ema_fixture(), ["EGFR"])
        self.assertEqual([row["therapy_name"] for row in rows], ["Old-EGFR"])


class WriterAndCliTests(unittest.TestCase):
    def test_csv_is_accepted_by_regulatory_provider(self):
        rows = parse_fda_html(FDA_FIXTURE, ["EGFR"])
        with tempfile.TemporaryDirectory() as directory:
            output = write_table(Path(directory) / "fda.csv", rows)
            records = RegulatoryTableTherapyProvider("fda", output).collect(
                TherapyContext(
                    gene="EGFR",
                    protein_changes=(),
                    disease_ids=(),
                    disease_names=(),
                )
            )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].therapy_status, "label_biomarker")
        self.assertEqual(records[0].molecular_profile, "EGFR (L858R)")

    def test_cli_uses_local_snapshots_without_network(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genes = root / "genes.txt"
            genes.write_text("EGFR\nMET\n", encoding="utf-8")
            fda = root / "fda.html"
            fda.write_text(FDA_FIXTURE, encoding="utf-8")
            ema = root / "ema.json"
            ema.write_text(json.dumps(ema_fixture()), encoding="utf-8")
            output = root / "out"
            code = main(
                [
                    "--gene-list",
                    str(genes),
                    "--fda-input",
                    str(fda),
                    "--ema-input",
                    str(ema),
                    "--output-dir",
                    str(output),
                ]
            )
            self.assertEqual(code, 0)
            self.assertTrue((output / "fda_biomarkers.csv").is_file())
            self.assertTrue((output / "ema_labels.csv").is_file())
            environment = (output / "regulatory_tables.env").read_text(
                encoding="utf-8"
            )
            self.assertIn("export FDA_FILE=", environment)
            self.assertIn("export EMA_FILE=", environment)
            with (output / "ema_labels.csv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)


if __name__ == "__main__":
    unittest.main()
