from __future__ import annotations

import unittest
from typing import Any

from rare_oncology.protein_annotations import UniProtProteinClient


class StubHttp:
    def __init__(self, payload: dict[str, Any]):
        self.payload = payload
        self.params: dict[str, object] | None = None

    def get_json(self, url: str, params: dict[str, Any]) -> dict[str, Any]:
        del url
        self.params = params
        return self.payload


UNIPROT_ACVR1 = {
    "results": [
        {
            "primaryAccession": "Q04771",
            "sequence": {"length": 509, "value": "A" * 509},
            "features": [
                {
                    "type": "Signal",
                    "location": {
                        "start": {"value": 1},
                        "end": {"value": 20},
                    },
                },
                {
                    "type": "Topological domain",
                    "description": "Extracellular",
                    "location": {
                        "start": {"value": 21},
                        "end": {"value": 180},
                    },
                },
                {
                    "type": "Transmembrane",
                    "description": "Helical",
                    "location": {
                        "start": {"value": 181},
                        "end": {"value": 203},
                    },
                },
                {
                    "type": "Domain",
                    "description": "Protein kinase",
                    "location": {
                        "start": {"value": 208},
                        "end": {"value": 502},
                    },
                },
                {
                    "type": "Helix",
                    "location": {
                        "start": {"value": 198},
                        "end": {"value": 205},
                    },
                },
                {
                    "type": "Beta strand",
                    "location": {
                        "start": {"value": 194},
                        "end": {"value": 196},
                    },
                },
                {
                    "type": "Binding site",
                    "location": {
                        "start": {"value": 300},
                        "end": {"value": 300},
                    },
                },
            ],
        }
    ]
}


class UniProtProteinClientTests(unittest.TestCase):
    def test_parses_domains_secondary_structure_and_exact_gene_query(self):
        http = StubHttp(UNIPROT_ACVR1)
        annotation = UniProtProteinClient(http=http).fetch("acvr1")

        self.assertIsNotNone(annotation)
        assert annotation is not None
        self.assertEqual(annotation.accession, "Q04771")
        self.assertEqual(annotation.length, 509)
        self.assertEqual(annotation.sequence, "A" * 509)
        self.assertEqual(
            [(feature.kind, feature.label) for feature in annotation.domains],
            [
                ("signal", "Signal peptide"),
                ("topology", "Extracellular"),
                ("transmembrane", "Transmembrane"),
                ("domain", "Protein kinase"),
            ],
        )
        self.assertEqual(
            [feature.kind for feature in annotation.secondary_structure],
            ["strand", "helix"],
        )
        params = http.params
        self.assertIsNotNone(params)
        assert params is not None
        self.assertIn("gene_exact:ACVR1", str(params["query"]))
        self.assertIn("organism_id:9606", str(params["query"]))
        self.assertIn("reviewed:true", str(params["query"]))

    def test_returns_none_when_uniprot_has_no_reviewed_entry(self):
        self.assertIsNone(UniProtProteinClient(http=StubHttp({"results": []})).fetch("X"))


if __name__ == "__main__":
    unittest.main()
