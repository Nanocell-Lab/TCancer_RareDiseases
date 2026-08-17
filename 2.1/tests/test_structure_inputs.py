from __future__ import annotations

import json
import sqlite3
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from rare_oncology.cli import main
from rare_oncology.models import AnalysisResult, SharedMutation
from rare_oncology.protein_annotations import ProteinAnnotation
from rare_oncology.structure_inputs import (
    StructureInputError,
    generate_variant_structure_inputs,
    read_protein_sequence,
)
from rare_oncology.structure_cli import main as structure_main


def shared(change: str) -> SharedMutation:
    return SharedMutation(
        gene="TESTR",
        protein_change=change,
        diseases=("Rare receptor disorder",),
        orphanet_codes=("1",),
        clinvar_ids=("1",),
        clinvar_accessions=("VCV1",),
        germline_classifications=("Pathogenic",),
        canonical_spdis=(),
        tcga_projects=("TCGA-TEST",),
        cancer_patient_count=1,
        cancer_event_count=1,
        genomic_match_count=0,
        genomic_confirmations=(),
        ssm_ids=("ssm-1",),
        match_level="candidate_protein_match",
    )


def result(*changes: str) -> AnalysisResult:
    mutations = tuple(shared(change) for change in changes)
    return AnalysisResult(
        gene="TESTR",
        generated_at="2026-08-16T00:00:00+00:00",
        gdc_program="TCGA",
        gdc_projects=(),
        min_cancer_patients=1,
        excluded_conditions=(),
        clinvar_record_limit=None,
        allowed_germline_classifications=("Pathogenic",),
        gdc_data_release=None,
        rare_variant_count=len(mutations),
        cancer_observation_count=len(mutations),
        shared_mutation_count=len(mutations),
        shared_mutations=mutations,
        mutation_plot_data=(),
    )


class StructureInputTests(unittest.TestCase):
    def test_writes_equivalent_af3_json_and_boltz2_yaml(self):
        with tempfile.TemporaryDirectory() as directory:
            generated = generate_variant_structure_inputs(
                result("L2R"),
                directory,
                protein_sequence="ALC",
                protein_accession="P12345",
                model_seeds=(7, 11),
            )
            af3 = next(
                artifact
                for artifact in generated.artifacts
                if artifact.platform == "alphafold3"
            )
            boltz = next(
                artifact
                for artifact in generated.artifacts
                if artifact.platform == "boltz2"
            )
            af3_payload = json.loads(Path(af3.path).read_text(encoding="utf-8"))
            boltz_text = Path(boltz.path).read_text(encoding="utf-8")
            manifest = json.loads(
                Path(generated.manifest_path).read_text(encoding="utf-8")
            )

        self.assertEqual(af3_payload["dialect"], "alphafold3")
        self.assertEqual(af3_payload["version"], 4)
        self.assertEqual(af3_payload["modelSeeds"], [7, 11])
        self.assertEqual(
            af3_payload["sequences"][0]["protein"]["sequence"], "ARC"
        )
        self.assertIn("version: 1", boltz_text)
        self.assertIn('sequence: "ARC"', boltz_text)
        self.assertIn("boltz predict INPUT.yaml --use_msa_server", boltz_text)
        self.assertEqual(af3.sequence_sha256, boltz.sequence_sha256)
        self.assertEqual(manifest["accession"], "P12345")
        self.assertEqual(len(manifest["artifacts"]), 2)

    def test_reference_mismatch_and_non_missense_are_skipped(self):
        with tempfile.TemporaryDirectory() as directory:
            generated = generate_variant_structure_inputs(
                result("A2R", "L2*", "L2_L3del"),
                directory,
                protein_sequence="ALC",
            )
            output_files = sorted(Path(directory).glob("*"))
        self.assertEqual(generated.artifacts, ())
        self.assertEqual(len(generated.skipped), 3)
        self.assertIn("UniProt contains L2", generated.skipped[0].reason)
        self.assertEqual(len(output_files), 1)
        self.assertTrue(output_files[0].name.endswith("_structure_inputs.json"))

    def test_explicit_variants_override_shared_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            generated = generate_variant_structure_inputs(
                result("L2R"),
                directory,
                variants=("p.C3A",),
                protein_sequence="ALC",
            )
        self.assertEqual(
            {item.protein_change for item in generated.artifacts}, {"C3A"}
        )

    def test_reads_single_fasta_and_rejects_noncanonical_residues(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "protein.fasta"
            fasta.write_text(">sp|P12345|TESTR\nAL\nC\n", encoding="utf-8")
            self.assertEqual(read_protein_sequence(fasta), "ALC")
            fasta.write_text(">bad\nAXC\n", encoding="utf-8")
            with self.assertRaisesRegex(
                StructureInputError, "non-canonical residues"
            ):
                read_protein_sequence(fasta)
            fasta.write_text(">one\nALC\n>two\nARC\n", encoding="utf-8")
            with self.assertRaisesRegex(StructureInputError, "a single sequence"):
                read_protein_sequence(fasta)

    def test_uniprot_sequence_is_reused_from_sqlite_in_offline_mode(self):
        class Client:
            def __init__(self):
                self.calls = 0

            def fetch(self, gene):
                self.calls += 1
                return ProteinAnnotation(
                    gene=gene,
                    accession="P12345",
                    length=3,
                    domains=(),
                    secondary_structure=(),
                    source_url="https://rest.uniprot.org/uniprotkb/P12345.json",
                    sequence="ALC",
                )

        class NoWeb:
            def fetch(self, gene: str) -> ProteinAnnotation | None:
                del gene
                raise AssertionError("no debe consultar UniProt")

        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "runs.sqlite3"
            client = Client()
            first = generate_variant_structure_inputs(
                result("L2R"),
                Path(directory) / "first",
                protein_client=client,
                sequence_database=database,
            )
            second = generate_variant_structure_inputs(
                result("C3A"),
                Path(directory) / "second",
                protein_client=NoWeb(),
                sequence_database=database,
                offline=True,
            )
            with sqlite3.connect(database) as connection:
                hits = connection.execute(
                    "SELECT cache_hit_count FROM protein_sequence_cache "
                    "WHERE gene = 'TESTR'"
                ).fetchone()[0]
        self.assertEqual(client.calls, 1)
        self.assertEqual(first.accession, "P12345")
        self.assertEqual(second.accession, "P12345")
        self.assertEqual(hits, 1)

    def test_offline_sequence_cache_miss_is_explicit(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(
                StructureInputError, "No UniProt sequence is cached"
            ):
                generate_variant_structure_inputs(
                    result("L2R"),
                    Path(directory) / "out",
                    sequence_database=Path(directory) / "runs.sqlite3",
                    offline=True,
                )

    def test_no_variants_or_invalid_seeds_fail_before_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(StructureInputError, "No shared or explicit variants"):
                generate_variant_structure_inputs(
                    result(), directory, protein_sequence="ALC"
                )
            with self.assertRaisesRegex(ValueError, "at least one seed"):
                generate_variant_structure_inputs(
                    result("L2R"),
                    directory,
                    protein_sequence="ALC",
                    model_seeds=(),
                )

    def test_cli_generates_both_platform_inputs_from_local_sequence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "protein.fasta"
            fasta.write_text(">TESTR\nALC\n", encoding="utf-8")
            output = root / "results"
            with patch(
                "rare_oncology.cli.find_shared_mutations",
                return_value=result("L2R"),
            ), redirect_stdout(StringIO()) as stdout:
                exit_code = main(
                    [
                        "TESTR",
                        "--no-database",
                        "--structure-variant",
                        "L2R",
                        "--protein-sequence-file",
                        str(fasta),
                        "--output-dir",
                        str(output),
                    ]
                )
            files = {path.name for path in (output / "structures").iterdir()}
        self.assertEqual(exit_code, 0)
        self.assertIn("testr_l2r_alphafold3.json", files)
        self.assertIn("testr_l2r_boltz2.yaml", files)
        self.assertIn("Structures: 2 files", stdout.getvalue())

    def test_focused_cli_uses_stored_report_and_local_fasta(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            report = root / "report.json"
            report.write_text(json.dumps(result("L2R").to_dict()), encoding="utf-8")
            fasta = root / "protein.fasta"
            fasta.write_text(">TESTR\nALC\n", encoding="utf-8")
            output = root / "focused"
            self.assertEqual(
                structure_main(
                    [
                        "--report-json",
                        str(report),
                        "--database",
                        str(root / "master.sqlite3"),
                        "--output-dir",
                        str(output),
                        "--variant",
                        "L2R",
                        "--protein-sequence-file",
                        str(fasta),
                        "--offline",
                    ]
                ),
                0,
            )
            files = {path.name for path in output.iterdir()}
        self.assertIn("testr_l2r_alphafold3.json", files)
        self.assertIn("testr_l2r_boltz2.yaml", files)
        self.assertIn("testr_structure_inputs.json", files)


if __name__ == "__main__":
    unittest.main()
