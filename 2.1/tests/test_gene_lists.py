from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ReceptorGeneListTests(unittest.TestCase):
    @staticmethod
    def _genes(filename: str) -> list[str]:
        return (ROOT / "gene_lists" / filename).read_text(
            encoding="utf-8"
        ).splitlines()

    def test_rtk_list_has_58_sorted_unique_symbols(self):
        genes = self._genes("receptor_tyrosine_kinases.txt")
        self.assertEqual(len(genes), 58)
        self.assertEqual(genes, sorted(set(genes)))
        self.assertTrue(
            {"EGFR", "ERBB2", "FGFR1", "MET", "RET", "ALK", "STYK1"}
            <= set(genes)
        )

    def test_kinase_linked_cytokine_receptors_are_complementary(self):
        genes = self._genes("kinase_linked_cytokine_receptors.txt")
        rtk = set(self._genes("receptor_tyrosine_kinases.txt"))
        tgfb = set(self._genes("tgfb_receptors.txt"))

        self.assertEqual(len(genes), 43)
        self.assertEqual(genes, sorted(set(genes)))
        self.assertFalse(set(genes) & (rtk | tgfb))
        self.assertTrue(
            {
                "CSF3R",
                "EPOR",
                "IFNAR1",
                "IL2RG",
                "IL6ST",
                "MPL",
                "PRLR",
            }
            <= set(genes)
        )
        self.assertFalse({"JAK1", "JAK2", "JAK3", "TYK2"} & set(genes))

    def test_maximum_runner_is_executable_and_valid_bash(self):
        script = ROOT / "scripts" / "run_rtk_maximum.sh"
        self.assertTrue(os.access(script, os.X_OK))
        completed = subprocess.run(
            ["bash", "-n", str(script)],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_runner_requests_all_data_layers_and_plot_formats(self):
        script = (ROOT / "scripts" / "run_rtk_maximum.sh").read_text(
            encoding="utf-8"
        )
        for required in (
            "--program ALL",
            "--classification ALL",
            "--source alphamissense",
            "--source gnomad",
            "--source panelapp",
            "--therapy-source civic",
            "--therapy-source oncokb",
            "--therapy-source clinicaltrials",
            "--therapy-source opentargets",
            "--plot-format png",
            "--plot-format svg",
            "--max-lollipop-residues",
            "--max-lollipop-labels",
            "--max-heatmap-mutations",
            "--max-rare-disease-mutations",
            "--string-network",
            "--string-partners",
            "--string-score",
            "--string-network-type",
        ):
            with self.subTest(required=required):
                self.assertIn(required, script)


if __name__ == "__main__":
    unittest.main()
