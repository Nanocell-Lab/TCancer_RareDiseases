from __future__ import annotations

import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest import mock

from rare_oncology.models import (
    AnalysisResult,
    MutationPlotRecord,
    SharedMutation,
)
from rare_oncology.plot_data import (
    heatmap_data,
    lollipop_label_residues,
    lollipop_records,
    lollipop_residues,
    rare_disease_heatmap_data,
    recurrence_records,
    venn_counts,
)
from rare_oncology.plots import (
    create_dashboard,
    create_heatmap_plot,
    create_lollipop_plot,
    create_rare_disease_heatmap_plot,
    create_recurrence_plot,
    create_venn_plot,
    generate_gene_plots,
)
from rare_oncology.http import ApiError
from rare_oncology.protein_annotations import ProteinAnnotation, ProteinFeature


def _shared(change: str, patients: int) -> SharedMutation:
    return SharedMutation(
        gene="ACVR1",
        protein_change=change,
        diseases=("Progressive myositis ossificans",),
        orphanet_codes=("337",),
        clinvar_ids=(f"clinvar-{change}",),
        clinvar_accessions=(f"VCV-{change}",),
        germline_classifications=("Pathogenic",),
        canonical_spdis=(),
        tcga_projects=("TCGA-UCEC",),
        cancer_patient_count=patients,
        cancer_event_count=1,
        genomic_match_count=1,
        genomic_confirmations=(f"VCV-{change}|ssm-{change}|chr2:1C>T",),
        ssm_ids=(f"ssm-{change}",),
        match_level="confirmed_genomic_and_protein",
    )


def _record(
    change: str,
    position: int,
    *,
    rare: bool,
    oncology: bool,
    patients: int = 0,
    projects: tuple[tuple[str, int], ...] = (),
) -> MutationPlotRecord:
    return MutationPlotRecord(
        gene="ACVR1",
        protein_change=change,
        amino_acid_position=position,
        in_rare_disease=rare,
        in_oncology=oncology,
        diseases=("Progressive myositis ossificans",) if rare else (),
        orphanet_codes=("337",) if rare else (),
        germline_classifications=("Pathogenic",) if rare else (),
        cancer_patient_count=patients,
        cancer_event_count=1 if oncology else 0,
        project_patient_counts=projects,
        match_level=("confirmed_genomic_and_protein" if rare and oncology else None),
    )


def _analysis_result() -> AnalysisResult:
    shared = (_shared("R206H", 7), _shared("G356D", 2))
    plot_data = (
        _record(
            "R206H",
            206,
            rare=True,
            oncology=True,
            patients=7,
            projects=(("TCGA-BRCA", 1), ("TCGA-UCEC", 6)),
        ),
        _record(
            "G356D",
            356,
            rare=True,
            oncology=True,
            patients=2,
            projects=(("TCGA-UCEC", 2),),
        ),
        _record("Q100*", 100, rare=True, oncology=False),
        _record(
            "K50N",
            50,
            rare=False,
            oncology=True,
            patients=1,
            projects=(("TCGA-LUAD", 1),),
        ),
        _record(
            "D200N",
            200,
            rare=False,
            oncology=True,
            patients=3,
            projects=(("TCGA-BRCA", 2), ("TCGA-LUAD", 1)),
        ),
    )
    return AnalysisResult(
        gene="ACVR1",
        generated_at="2026-08-11T00:00:00+00:00",
        gdc_program="TCGA",
        gdc_projects=(),
        min_cancer_patients=1,
        excluded_conditions=(),
        clinvar_record_limit=None,
        allowed_germline_classifications=("Pathogenic",),
        gdc_data_release="Data Release fixture",
        rare_variant_count=3,
        cancer_observation_count=13,
        shared_mutation_count=2,
        shared_mutations=shared,
        mutation_plot_data=plot_data,
    )


def _protein_annotation() -> ProteinAnnotation:
    return ProteinAnnotation(
        gene="ACVR1",
        accession="Q04771",
        length=509,
        domains=(
            ProteinFeature("domain", 178, 207, "GS"),
            ProteinFeature("domain", 208, 502, "Protein kinase"),
        ),
        secondary_structure=(
            ProteinFeature("helix", 180, 184, "Helix"),
            ProteinFeature("strand", 194, 196, "β-strand"),
            ProteinFeature("helix", 198, 205, "Helix"),
        ),
        source_url="https://rest.uniprot.org/uniprotkb/Q04771.json",
    )


class PlotDataTests(unittest.TestCase):
    def setUp(self):
        self.result = _analysis_result()

    def test_venn_uses_unique_normalized_mutations(self):
        counts = venn_counts(self.result)
        self.assertEqual((counts.rare_only, counts.shared, counts.oncology_only), (1, 2, 2))

    def test_lollipop_records_are_positioned_and_shared_first_before_limit(self):
        records = lollipop_records(self.result, max_mutations=2)
        self.assertEqual({record.protein_change for record in records}, {"R206H", "G356D"})

    def test_lollipop_prioritizes_rare_over_oncology_only_after_shared(self):
        records = lollipop_records(self.result, max_mutations=3)
        self.assertEqual(
            {record.protein_change for record in records},
            {"R206H", "G356D", "Q100*"},
        )

    def test_residue_overview_is_balanced_and_excludes_synonymous(self):
        dense = tuple(
            _record(
                f"A{position}V",
                position,
                rare=True,
                oncology=True,
                patients=1,
            )
            for position in range(100, 220)
        ) + tuple(
            _record(f"G{position}D", position, rare=True, oncology=False)
            for position in range(300, 500)
        ) + (
            _record("E746_A750del", 746, rare=False, oncology=True, patients=115),
            _record("L858R", 858, rare=False, oncology=True, patients=182),
            _record("R2=", 2, rare=True, oncology=False),
        )
        result = replace(self.result, mutation_plot_data=dense)

        residues = lollipop_residues(result, max_residues=80, max_position=1210)

        changes = {
            change for residue in residues for change in residue.protein_changes
        }
        self.assertIn("L858R", changes)
        self.assertIn("E746_A750del", changes)
        self.assertNotIn("R2=", changes)
        self.assertTrue(any(residue.category == "shared" for residue in residues))
        self.assertTrue(any(residue.category == "oncology" for residue in residues))
        self.assertTrue(any(residue.category == "rare" for residue in residues))

    def test_residue_overview_aggregates_changes_at_one_position(self):
        result = replace(
            self.result,
            mutation_plot_data=(
                _record("R206H", 206, rare=True, oncology=True, patients=7),
                _record("R206C", 206, rare=True, oncology=False),
                _record("R206S", 206, rare=False, oncology=True, patients=2),
            ),
        )
        residues = lollipop_residues(result)
        self.assertEqual(len(residues), 1)
        self.assertEqual(residues[0].position, 206)
        self.assertEqual(
            (residues[0].rare_count, residues[0].oncology_count, residues[0].shared_count),
            (1, 1, 1),
        )
        self.assertEqual(residues[0].cancer_patient_count, 7)

    def test_lollipop_labels_are_hard_limited(self):
        residues = lollipop_residues(self.result, max_residues=5)
        labels = lollipop_label_residues(residues, max_labels=3)
        self.assertEqual(len(labels), 3)

    def test_recurrence_is_sorted_by_patient_count(self):
        records = recurrence_records(self.result)
        self.assertEqual([record.protein_change for record in records], ["R206H", "G356D"])

    def test_heatmap_counts_projects(self):
        data = heatmap_data(self.result, max_mutations=5)
        self.assertEqual(data.projects, ("TCGA-BRCA", "TCGA-LUAD", "TCGA-UCEC"))
        row = data.values[data.mutations.index("R206H")]
        self.assertEqual(row, (1, 0, 6))

    def test_rare_disease_heatmap_prioritizes_shared_within_cap(self):
        shared = replace(
            _record("R206H", 206, rare=True, oncology=True, patients=1),
            diseases=("Disease A",),
            disease_variant_counts=(("Disease A", 1),),
        )
        abundant_rare = replace(
            _record("Q100*", 100, rare=True, oncology=False),
            diseases=("Disease A", "Disease B"),
            disease_variant_counts=(("Disease A", 8), ("Disease B", 4)),
        )
        other_rare = replace(
            _record("A50V", 50, rare=True, oncology=False),
            diseases=("Disease C",),
            disease_variant_counts=(("Disease C", 6),),
        )
        result = replace(
            self.result,
            mutation_plot_data=(shared, abundant_rare, other_rare),
        )

        data = rare_disease_heatmap_data(result, max_mutations=2)

        self.assertEqual(data.mutations, ("R206H", "Q100*"))
        self.assertEqual(data.shared_mutations, (True, False))
        self.assertEqual(data.diseases, ("Disease A", "Disease B"))
        self.assertEqual(data.values, ((1, 8), (0, 4)))

    def test_rare_disease_heatmap_has_hard_limit_of_thirty(self):
        with self.assertRaisesRegex(ValueError, "between 1 and 30"):
            rare_disease_heatmap_data(self.result, max_mutations=31)


class PlotRenderingTests(unittest.TestCase):
    def setUp(self):
        self.result = _analysis_result()
        self.annotation = _protein_annotation()
        self.figures = []

    def tearDown(self):
        import matplotlib.pyplot as plt

        for figure in self.figures:
            plt.close(figure)

    def test_plot_semantics(self):
        lollipop = create_lollipop_plot(
            self.result,
            protein_annotation=self.annotation,
        )
        venn = create_venn_plot(self.result)
        recurrence = create_recurrence_plot(self.result)
        heatmap = create_heatmap_plot(self.result)
        rare_heatmap = create_rare_disease_heatmap_plot(self.result)
        dashboard = create_dashboard(
            self.result,
            protein_annotation=self.annotation,
        )
        self.figures.extend(
            [lollipop, venn, recurrence, heatmap, rare_heatmap, dashboard]
        )

        self.assertIn("variant protein map", lollipop.axes[0].get_title())
        self.assertIn("UniProt Q04771", lollipop.axes[0].get_title())
        track_labels = {text.get_text() for text in lollipop.axes[0].texts}
        self.assertTrue({"Variants", "Domains", "Sec. structure"} <= track_labels)
        self.assertTrue(any(label.startswith("p.R206H") for label in track_labels))
        self.assertEqual(
            [tick.get_text() for tick in lollipop.axes[0].get_xticklabels()][::5],
            ["1", "509"],
        )
        self.assertEqual(len(venn.axes[0].patches), 2)
        self.assertEqual(
            [patch.get_height() for patch in recurrence.axes[0].patches],
            [7, 2],
        )
        self.assertEqual(len(heatmap.axes[0].images), 1)
        self.assertEqual(len(rare_heatmap.axes[0].images), 1)
        self.assertTrue(
            any(
                tick.get_text().startswith("[C] p.")
                for tick in rare_heatmap.axes[0].get_xticklabels()
            )
        )
        self.assertGreaterEqual(len(dashboard.axes), 4)

    def test_generate_gene_plots_writes_nonempty_png_and_svg(self):
        with tempfile.TemporaryDirectory() as directory:
            written = generate_gene_plots(
                self.result,
                directory,
                protein_annotation=self.annotation,
            )
            self.assertEqual(
                set(written),
                {
                    "lollipop",
                    "venn",
                    "recurrence",
                    "heatmap",
                    "rare_disease_heatmap",
                    "dashboard",
                },
            )
            paths = [path for group in written.values() for path in group]
            self.assertEqual(len(paths), 12)
            for path in paths:
                self.assertTrue(path.is_file())
                self.assertGreater(path.stat().st_size, 100)
            png = Path(directory) / "acvr1_lollipop.png"
            self.assertEqual(png.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
            svg = Path(directory) / "acvr1_lollipop.svg"
            self.assertIn("<svg", svg.read_text(encoding="utf-8")[:1000])

    def test_invalid_limit_creates_no_directory_and_leaks_no_figure(self):
        import matplotlib.pyplot as plt

        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "plots"
            with self.assertRaises(ValueError):
                generate_gene_plots(
                    self.result,
                    destination,
                    max_lollipop_mutations=0,
                )
            self.assertFalse(destination.exists())
            self.assertEqual(plt.get_fignums(), [])

    def test_public_renderer_closes_figure_when_drawing_fails(self):
        import matplotlib.pyplot as plt

        with mock.patch(
            "rare_oncology.plots._draw_lollipop",
            side_effect=RuntimeError("fixture"),
        ):
            with self.assertRaises(RuntimeError):
                create_lollipop_plot(
                    self.result,
                    protein_annotation=self.annotation,
                )
        self.assertEqual(plt.get_fignums(), [])

    def test_variant_outside_uniprot_length_is_omitted_and_reported(self):
        outside = _record("A600V", 600, rare=True, oncology=False)
        result = replace(
            self.result,
            mutation_plot_data=self.result.mutation_plot_data + (outside,),
        )
        figure = create_lollipop_plot(
            result,
            protein_annotation=self.annotation,
        )
        self.figures.append(figure)

        labels = {text.get_text() for text in figure.axes[0].texts}
        self.assertNotIn("p.A600V", labels)
        self.assertIn("1 out-of-range change omitted", figure.axes[0].get_title())
        self.assertEqual(figure.axes[0].get_xlim()[1], 510)

    def test_dense_lollipop_keeps_exact_bounds_and_oncology_hotspots(self):
        dense = tuple(
            _record(
                f"A{position}V",
                position,
                rare=True,
                oncology=True,
                patients=1,
            )
            for position in range(100, 220)
        ) + tuple(
            _record(f"G{position}D", position, rare=True, oncology=False)
            for position in range(300, 700)
        ) + (
            _record("E746_A750del", 746, rare=False, oncology=True, patients=115),
            _record("L858R", 858, rare=False, oncology=True, patients=182),
        )
        result = replace(self.result, gene="EGFR", mutation_plot_data=dense)
        annotation = ProteinAnnotation(
            gene="EGFR",
            accession="P00533",
            length=1210,
            domains=(ProteinFeature("domain", 712, 979, "Protein kinase"),),
            secondary_structure=(),
            source_url="https://rest.uniprot.org/uniprotkb/P00533.json",
        )
        figure = create_lollipop_plot(
            result,
            max_mutations=80,
            max_labels=15,
            protein_annotation=annotation,
        )
        self.figures.append(figure)
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        tight = figure.get_tightbbox(renderer)
        labels = [
            text.get_text()
            for text in figure.axes[0].texts
            if text.get_text().startswith("p.")
        ]

        self.assertLessEqual(tight.width, figure.get_figwidth() * 1.15)
        self.assertLessEqual(len(labels), 15)
        self.assertTrue(any(label.startswith("p.L858R") for label in labels))
        self.assertTrue(any(label.startswith("p.E746_A750del") for label in labels))
        self.assertEqual(tuple(round(value) for value in figure.axes[0].get_xlim()), (1, 1211))

    def test_uniprot_api_error_uses_explicit_inferred_fallback(self):
        with mock.patch(
            "rare_oncology.plots.fetch_protein_annotation",
            side_effect=ApiError("fixture"),
        ):
            with self.assertWarnsRegex(RuntimeWarning, "Fallo consultando|fixture"):
                figure = create_lollipop_plot(self.result)
        self.figures.append(figure)
        self.assertIn("inferred length (no UniProt)", figure.axes[0].get_title())

    def test_provided_offline_annotation_never_queries_uniprot(self):
        with mock.patch("rare_oncology.plots.fetch_protein_annotation") as fetch:
            figure = create_lollipop_plot(
                self.result,
                protein_annotation=ProteinAnnotation(
                    gene=self.result.gene,
                    accession=None,
                    length=900,
                    domains=(),
                    secondary_structure=(),
                    source_url=None,
                ),
            )
        self.figures.append(figure)
        fetch.assert_not_called()

    def test_missing_uniprot_entry_warns_and_labels_fallback(self):
        with mock.patch(
            "rare_oncology.plots.fetch_protein_annotation",
            return_value=None,
        ):
            with self.assertWarnsRegex(RuntimeWarning, "No reviewed human UniProt entry"):
                figure = create_lollipop_plot(self.result)
        self.figures.append(figure)
        self.assertIn("inferred length (no UniProt)", figure.axes[0].get_title())


if __name__ == "__main__":
    unittest.main()
