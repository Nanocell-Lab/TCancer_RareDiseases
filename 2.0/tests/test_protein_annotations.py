from __future__ import annotations

import unittest

from rare_oncology.protein_annotations import UniProtProteinClient


class StubHttp:
    def __init__(self, payload):
        self.payload = payload
        self.params = None

    def get_json(self, _url, params):
        self.params = params
        return self.payload


UNIPROT_ACVR1 = {
    "results": [
        {
            "primaryAccession": "Q04771",
            "sequence": {"length": 509},
            "features": [
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
        self.assertEqual(annotation.domains[0].label, "Protein kinase")
        self.assertEqual(
            [feature.kind for feature in annotation.secondary_structure],
            ["strand", "helix"],
        )
        self.assertIn("gene_exact:ACVR1", http.params["query"])
        self.assertIn("organism_id:9606", http.params["query"])
        self.assertIn("reviewed:true", http.params["query"])

    def test_returns_none_when_uniprot_has_no_reviewed_entry(self):
        self.assertIsNone(UniProtProteinClient(http=StubHttp({"results": []})).fetch("X"))


if __name__ == "__main__":
    unittest.main()
