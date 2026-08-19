from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from rare_oncology_dashboard.runner import (
    AnalysisRunner,
    RunOptions,
    RunnerError,
    StringRunOptions,
    StructureRunOptions,
    normalize_gene,
)


class RunnerTests(unittest.TestCase):
    def test_command_links_to_21_and_never_contains_token(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / "2.1" / ".venv" / "bin" / "rare-oncology"
            executable.parent.mkdir(parents=True)
            executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
            master = root / "2.1" / "results" / "dashboard_master.sqlite3"
            runner = AnalysisRunner(
                root,
                master,
                environment={"ONCOKB_TOKEN": "TOP-SECRET"},
            )
            command = runner.build_command(
                "egfr",
                RunOptions(include_oncokb=True),
            )
            self.assertEqual(command[1], "EGFR")
            self.assertIn(str(master.resolve()), command)
            self.assertIn("oncokb", command)
            self.assertNotIn("TOP-SECRET", " ".join(command))
            self.assertNotIn("--plots", command)
            self.assertIn("--string-network", command)
            self.assertEqual(command[command.index("--string-score") + 1], "700")

    def test_oncokb_requires_environment_token(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / "2.1" / ".venv" / "bin" / "rare-oncology"
            executable.parent.mkdir(parents=True)
            executable.touch()
            runner = AnalysisRunner(root, root / "master.sqlite3", environment={})
            with self.assertRaises(RunnerError):
                runner.build_command("EGFR", RunOptions(include_oncokb=True))

    def test_clone_tables_are_added_to_therapy_runs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / "2.1" / ".venv" / "bin" / "rare-oncology"
            executable.parent.mkdir(parents=True)
            executable.touch()
            regulatory = root / "2.1" / "bootstrap" / "regulatory"
            regulatory.mkdir(parents=True)
            fda = regulatory / "fda_biomarkers.csv"
            ema = regulatory / "ema_labels.csv"
            fda.write_text("gene\nEGFR\n", encoding="utf-8")
            ema.write_text("gene\nEGFR\n", encoding="utf-8")
            runner = AnalysisRunner(root, root / "master.sqlite3", environment={})

            command = runner.build_command("EGFR", RunOptions())

            self.assertEqual(command.count("--regulatory-file"), 2)
            self.assertIn(f"fda={fda.resolve()}", command)
            self.assertIn(f"ema={ema.resolve()}", command)

    def test_rejects_shell_metacharacters_in_gene(self) -> None:
        with self.assertRaises(ValueError):
            normalize_gene("EGFR; rm -rf /")

    def test_executes_fake_21_cli_without_a_shell(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / "2.1" / ".venv" / "bin" / "rare-oncology"
            executable.parent.mkdir(parents=True)
            executable.write_text("#!/bin/sh\nprintf 'fake run: %s\\n' \"$1\"\n", encoding="utf-8")
            executable.chmod(0o755)
            runner = AnalysisRunner(root, root / "2.1" / "results" / "master.sqlite3")
            outcome = runner.run(
                "TGFBR1",
                RunOptions(
                    include_therapies=False,
                    evidence_sources=(),
                    timeout_seconds=5,
                ),
            )
            self.assertEqual(outcome.returncode, 0)
            self.assertIn("fake run: TGFBR1", outcome.stdout)
            self.assertTrue(outcome.output_directory.is_dir())

    def test_builds_focused_string_command_without_full_analysis(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            python = root / "2.1" / ".venv" / "bin" / "python"
            python.parent.mkdir(parents=True)
            python.touch()
            master = root / "2.1" / "results" / "master.sqlite3"
            runner = AnalysisRunner(root, master, environment={"SECRET": "hidden"})
            command = runner.build_string_command(
                "tgfbr1",
                StringRunOptions(
                    required_score=800,
                    added_nodes=30,
                    network_type="physical",
                ),
            )
            self.assertEqual(
                command[:4],
                (str(python.resolve()), "-m", "rare_oncology.string_cli", "TGFBR1"),
            )
            self.assertEqual(command[command.index("--string-score") + 1], "800")
            self.assertEqual(command[command.index("--string-partners") + 1], "30")
            self.assertNotIn("rare-oncology", command)
            self.assertNotIn("hidden", " ".join(command))

    def test_builds_focused_structure_command_with_variants_and_seeds(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            python = root / "2.1" / ".venv" / "bin" / "python"
            python.parent.mkdir(parents=True)
            python.touch()
            runner = AnalysisRunner(root, root / "master.sqlite3", environment={})
            report = root / "report.json"
            sequence = root / "sequence.fasta"
            command = runner.build_structure_command(
                "egfr",
                StructureRunOptions(
                    variants=("L858R", "T790M"),
                    model_seeds=(1, 7),
                    offline=True,
                ),
                report_path=report,
                sequence_path=sequence,
            )
            self.assertEqual(command[2:4], ("rare_oncology.structure_cli", "--report-json"))
            self.assertEqual(command.count("--variant"), 2)
            self.assertEqual(command.count("--model-seed"), 2)
            self.assertIn("--offline", command)
            self.assertNotIn("--program", command)


if __name__ == "__main__":
    unittest.main()
