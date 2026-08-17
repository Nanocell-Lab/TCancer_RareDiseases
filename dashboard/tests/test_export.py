from __future__ import annotations

import unittest
import zipfile
from io import BytesIO

from rare_oncology_dashboard.export import files_to_zip, records_to_csv


class ExportTests(unittest.TestCase):
    def test_csv_neutralizes_formula_cells(self) -> None:
        payload = records_to_csv(
            [{"safe": "EGFR", "formula": "=CMD()", "hidden": "\t@SUM(A1:A2)"}]
        ).decode("utf-8-sig")
        self.assertIn("'=CMD()", payload)
        self.assertIn("'\t@SUM", payload)

    def test_structure_zip_is_flat_and_preserves_content(self) -> None:
        payload = files_to_zip([("path/input.json", b"{}"), ("input.yaml", b"v: 1")])
        with zipfile.ZipFile(BytesIO(payload)) as archive:
            self.assertEqual(set(archive.namelist()), {"input.json", "input.yaml"})
            self.assertEqual(archive.read("input.json"), b"{}")


if __name__ == "__main__":
    unittest.main()
