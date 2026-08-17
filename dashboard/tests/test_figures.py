from __future__ import annotations

import unittest

from rare_oncology_dashboard.figures import (
    heatmap_figure,
    lollipop_figure,
    string_evidence_heatmap_figure,
    string_network_figure,
    string_partner_scores_figure,
    residue_displacement_figure,
    structure_overlay_figure,
)


class FigureTests(unittest.TestCase):
    def test_lollipop_contains_baseline_stems_and_category(self) -> None:
        figure = lollipop_figure(
            [
                {
                    "position": 206,
                    "height": 3.0,
                    "category": "Shared",
                    "protein_changes": "R206H",
                    "mutation_count": 1,
                    "cancer_patient_count": 4,
                    "diseases": "FOP",
                    "projects": "TCGA-BRCA",
                    "confirmed": True,
                }
            ],
            "ACVR1",
        )
        payload = figure.to_dict()
        self.assertEqual(len(payload["data"]), 3)
        self.assertEqual(payload["data"][-1]["name"], "Shared")
        self.assertEqual(payload["data"][-1]["x"], [206])
        self.assertIn("ACVR1", str(figure.layout.title.text))
        figure.to_json()

    def test_rare_heatmap_marks_shared_columns(self) -> None:
        figure = heatmap_figure(
            {
                "diseases": ["Rare disease"],
                "mutations": ["R206H", "G356D"],
                "values": [[3, 1]],
                "shared": [True, False],
            },
            row_key="diseases",
            column_key="mutations",
            title="Rare",
            colorbar_title="Associations",
        )
        payload = figure.to_dict()
        self.assertEqual(payload["data"][0]["x"], ["★ p.R206H", "p.G356D"])
        self.assertEqual(payload["data"][0]["z"], [[3, 1]])

    def test_string_figures_expose_query_partner_and_scores(self) -> None:
        graph = {
            "query_id": "Q",
            "nodes": [
                {"string_id": "Q", "preferred_name": "EGFR", "is_query": True, "degree": 1},
                {"string_id": "P", "preferred_name": "GRB2", "is_query": False, "degree": 1},
            ],
            "edges": [{"string_id_a": "Q", "string_id_b": "P", "score": 0.9}],
        }
        network = string_network_figure(graph, "EGFR").to_dict()
        self.assertEqual(network["data"][-1]["name"], "EGFR")
        self.assertIn("GRB2", network["data"][-2]["text"])

        rows = [{
            "partner": "GRB2", "combined_score": 0.9,
            "experimental_score": 0.8, "database_score": 0.7,
            "coexpression_score": 0.6, "textmining_score": 0.5,
            "neighborhood_score": 0.0, "fusion_score": 0.0,
            "phylogenetic_score": 0.0,
        }]
        bars = string_partner_scores_figure(rows).to_dict()
        self.assertEqual(bars["data"][0]["y"], ["GRB2"])
        heatmap = string_evidence_heatmap_figure(rows).to_dict()
        self.assertEqual(heatmap["data"][0]["z"][0][0], 0.8)

    def test_structure_figures_show_aligned_models_and_residue_shift(self) -> None:
        comparison = {
            "wildtype_coordinates": [[0, 0, 0], [1, 0, 0]],
            "aligned_variant_coordinates": [[0, 0, 0], [1, 0.2, 0]],
            "residues": [
                {
                    "chain": "A", "residue_number": 1, "insertion_code": "",
                    "wildtype_residue": "ALA", "variant_residue": "ALA",
                    "residue_changed": False, "ca_displacement_angstrom": 0.0,
                },
                {
                    "chain": "A", "residue_number": 2, "insertion_code": "",
                    "wildtype_residue": "GLY", "variant_residue": "VAL",
                    "residue_changed": True, "ca_displacement_angstrom": 0.2,
                },
            ],
        }
        overlay = structure_overlay_figure(comparison, "EGFR").to_dict()
        self.assertEqual([trace["name"] for trace in overlay["data"]], ["Wild type", "Variant (aligned)"])
        shifts = residue_displacement_figure(comparison).to_dict()
        self.assertEqual(shifts["data"][0]["y"], [0.0, 0.2])


if __name__ == "__main__":
    unittest.main()
