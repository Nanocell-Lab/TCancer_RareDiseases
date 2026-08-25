from __future__ import annotations

import unittest

from rare_oncology.chemical_metadata import ChEMBLChemicalMetadataClient


class StubHttp:
    def __init__(self, payloads):
        self.payloads = iter(payloads)
        self.calls = []

    def get_json(self, url, params):
        self.calls.append((url, params))
        return next(self.payloads)


class ChEMBLChemicalMetadataTests(unittest.TestCase):
    def test_resolves_canonical_smiles_in_batches_and_caches_results(self):
        http = StubHttp(
            [
                {
                    "molecules": [
                        {
                            "molecule_chembl_id": "CHEMBL1",
                            "molecule_type": "Small molecule",
                            "molecule_structures": {
                                "canonical_smiles": "CCO"
                            },
                        }
                    ]
                },
                {
                    "molecules": [
                        {
                            "molecule_chembl_id": "CHEMBL2",
                            "molecule_type": "Antibody",
                            "molecule_structures": None,
                        }
                    ]
                },
            ]
        )
        client = ChEMBLChemicalMetadataClient(http=http, batch_size=1)

        first = client.resolve(("chembl2", "not-a-chembl-id", "CHEMBL1"))
        second = client.resolve(("CHEMBL1", "CHEMBL2"))

        self.assertEqual(list(first), ["CHEMBL1", "CHEMBL2"])
        self.assertEqual(first["CHEMBL1"].molecule_type, "Small molecule")
        self.assertEqual(first["CHEMBL1"].canonical_smiles, "CCO")
        self.assertIsNone(first["CHEMBL2"].canonical_smiles)
        self.assertEqual(second, first)
        self.assertEqual(len(http.calls), 2)
        self.assertEqual(
            http.calls[0][1],
            {"molecule_chembl_id__in": "CHEMBL1", "limit": 1},
        )

    def test_missing_molecule_is_negative_cached(self):
        http = StubHttp([{"molecules": []}])
        client = ChEMBLChemicalMetadataClient(http=http)

        self.assertEqual(client.resolve(("CHEMBL999",)), {})
        self.assertEqual(client.resolve(("CHEMBL999",)), {})
        self.assertEqual(len(http.calls), 1)


if __name__ == "__main__":
    unittest.main()
