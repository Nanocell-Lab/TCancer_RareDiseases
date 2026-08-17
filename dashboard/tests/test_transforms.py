from __future__ import annotations

import unittest

from rare_oncology_dashboard.transforms import (
    alphamissense_rows,
    lollipop_rows,
    oncology_heatmap,
    overlap_counts,
    rare_disease_heatmap,
    string_graph_data,
    string_partner_rows,
    structure_variant_rows,
    therapy_rows,
)


RESULT = {
    "mutation_plot_data": [
        {
            "protein_change": "R206H",
            "amino_acid_position": 206,
            "in_rare_disease": True,
            "in_oncology": True,
            "diseases": ["Fibrodysplasia ossificans progressiva"],
            "disease_variant_counts": [["Fibrodysplasia ossificans progressiva", 3]],
            "cancer_patient_count": 4,
            "project_patient_counts": [["TCGA-BRCA", 3], ["TCGA-COAD", 1]],
            "match_level": "confirmed_genomic_and_protein",
        },
        {
            "protein_change": "G356D",
            "amino_acid_position": 356,
            "in_rare_disease": True,
            "in_oncology": False,
            "diseases": ["Rare syndrome"],
            "disease_variant_counts": [["Rare syndrome", 2]],
            "cancer_patient_count": 0,
            "project_patient_counts": [],
            "match_level": None,
        },
        {
            "protein_change": "K50N",
            "amino_acid_position": 50,
            "in_rare_disease": False,
            "in_oncology": True,
            "diseases": [],
            "cancer_patient_count": 2,
            "project_patient_counts": [["TCGA-LUAD", 2]],
            "match_level": None,
        },
    ],
    "shared_mutations": [{"protein_change": "R206H"}],
    "evidence_records": [
        {
            "source": "alphamissense",
            "protein_change": "R206H",
            "prediction_score": 0.98,
            "prediction_class": "likely_pathogenic",
            "source_record_id": "am-1",
            "source_version": "v1",
        }
    ],
}

STRING_RESULT = {
    "mapped_string_id": "9606.ENSP_QUERY",
    "nodes": [
        {"string_id": "9606.ENSP_QUERY", "preferred_name": "EGFR", "is_query": True},
        {"string_id": "9606.ENSP_A", "preferred_name": "GRB2", "is_query": False},
        {"string_id": "9606.ENSP_B", "preferred_name": "ERBB2", "is_query": False},
    ],
    "edges": [
        {
            "string_id_a": "9606.ENSP_QUERY", "string_id_b": "9606.ENSP_A",
            "preferred_name_a": "EGFR", "preferred_name_b": "GRB2", "score": 0.96,
            "experimental_score": 0.8, "database_score": 0.7,
        },
        {
            "string_id_a": "9606.ENSP_B", "string_id_b": "9606.ENSP_QUERY",
            "preferred_name_a": "ERBB2", "preferred_name_b": "EGFR", "score": 0.91,
            "coexpression_score": 0.6, "textmining_score": 0.5,
        },
        {
            "string_id_a": "9606.ENSP_A", "string_id_b": "9606.ENSP_B",
            "preferred_name_a": "GRB2", "preferred_name_b": "ERBB2", "score": 0.88,
        },
    ],
}


class TransformTests(unittest.TestCase):
    def test_overlap_and_lollipop_preserve_three_categories(self) -> None:
        counts = {row["category"]: row["count"] for row in overlap_counts(RESULT)}
        self.assertEqual(counts, {"Rare disease only": 1, "Shared": 1, "Oncology only": 1})
        rows = lollipop_rows(RESULT)
        self.assertEqual([row["position"] for row in rows], [50, 206, 356])
        shared = next(row for row in rows if row["position"] == 206)
        self.assertTrue(shared["confirmed"])
        self.assertEqual(shared["cancer_patient_count"], 4)

    def test_oncology_heatmap_counts_projects(self) -> None:
        matrix = oncology_heatmap(RESULT)
        self.assertEqual(matrix["projects"], ["TCGA-BRCA", "TCGA-COAD", "TCGA-LUAD"])
        r206 = matrix["values"][matrix["mutations"].index("R206H")]
        self.assertEqual(r206, [3, 1, 0])

    def test_rare_heatmap_prioritizes_shared_mutations(self) -> None:
        matrix = rare_disease_heatmap(RESULT, top=1)
        self.assertEqual(matrix["mutations"], ["R206H"])
        self.assertEqual(matrix["shared"], [True])
        self.assertEqual(matrix["values"], [[3]])

    def test_alphamissense_marks_shared(self) -> None:
        rows = alphamissense_rows(RESULT)
        self.assertEqual(rows[0]["position"], 206)
        self.assertTrue(rows[0]["shared_with_oncology"])
        self.assertEqual(rows[0]["score"], 0.98)

    def test_therapy_layers_remain_separate(self) -> None:
        rows = therapy_rows(
            {
                "variant_therapy_evidence": [{"therapy_name": "A", "source": "CIViC"}],
                "regulatory_therapy": [{"therapy_name": "B", "source": "FDA"}],
                "trial_and_target_discovery": [{"therapy_name": "C", "source": "OT"}],
            }
        )
        self.assertEqual({row["layer"] for row in rows}, {
            "Variant-specific evidence", "Regulatory therapy", "Trial and target discovery"
        })

    def test_string_partners_are_direct_ranked_and_channel_safe(self) -> None:
        rows = string_partner_rows(STRING_RESULT)
        self.assertEqual([row["partner"] for row in rows], ["GRB2", "ERBB2"])
        self.assertEqual(rows[0]["combined_score"], 0.96)
        self.assertEqual(rows[0]["fusion_score"], 0.0)

    def test_string_graph_retains_direct_and_bounds_secondary_edges(self) -> None:
        graph = string_graph_data(STRING_RESULT, max_secondary_edges=0)
        self.assertEqual(graph["total_edges"], 3)
        self.assertEqual(graph["displayed_edges"], 2)
        self.assertEqual(len(graph["nodes"]), 3)
        query = next(node for node in graph["nodes"] if node["is_query"])
        self.assertEqual(query["degree"], 2)

    def test_structure_variants_include_only_simple_missense_and_prioritize_shared(self) -> None:
        payload = dict(RESULT)
        payload["mutation_plot_data"] = [
            *RESULT["mutation_plot_data"],
            {
                "protein_change": "E746_A750del",
                "in_rare_disease": True,
                "in_oncology": True,
            },
        ]
        rows = structure_variant_rows(payload)
        self.assertEqual(rows[0]["protein_change"], "R206H")
        self.assertNotIn("E746_A750del", {row["protein_change"] for row in rows})


if __name__ == "__main__":
    unittest.main()
