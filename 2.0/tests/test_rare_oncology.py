from __future__ import annotations

import unittest
import urllib.error
from unittest import mock

from rare_oncology.analysis import _aggregate_matches, find_shared_mutations
from rare_oncology.clients import GdcClient, NcbiClinVarClient
from rare_oncology.http import ApiError, JsonHttpClient
from rare_oncology.models import CancerObservation, RareVariant
from rare_oncology.normalization import normalize_protein_change
from rare_oncology.plot_data import venn_counts


class StubHttp:
    def __init__(self, *, get_responses=(), post_responses=()):
        self.get_responses = iter(get_responses)
        self.post_responses = iter(post_responses)

    def get_json(self, _url, _params):
        return next(self.get_responses)

    def post_json(self, _url, _payload):
        return next(self.post_responses)


CLINVAR_SEARCH = {
    "esearchresult": {"count": "2", "querykey": "1", "webenv": "fixture"}
}

GDC_STATUS = {"data_release": "Data Release fixture"}

CLINVAR_SUMMARY = {
    "result": {
        "uids": ["101", "102"],
        "101": {
            "uid": "101",
            "accession_version": "VCV000000101.1",
            "title": "NM_005228.5(EGFR):c.2573T>G (p.Leu858Arg)",
            "genes": [{"symbol": "EGFR"}],
            "protein_change": "L858R",
            "variation_set": [{"canonical_spdi": "NC_000007.14:55191821:T:G"}],
            "molecular_consequence_list": ["missense variant"],
            "germline_classification": {
                "description": "Pathogenic",
                "trait_set": [
                    {
                        "trait_name": "EGFR-associated rare disorder",
                        "trait_xrefs": [
                            {"db_source": "Orphanet", "db_id": "999999"}
                        ],
                    }
                ],
            },
        },
        "102": {
            "uid": "102",
            "accession_version": "VCV000000102.1",
            "title": "NM_005228.5(EGFR):c.3005T>A (p.Met1002Lys)",
            "genes": [{"symbol": "EGFR"}],
            "protein_change": "M1002K",
            "germline_classification": {
                "description": "Uncertain significance",
                "trait_set": [
                    {
                        "trait_name": "Condition without Orphanet mapping",
                        "trait_xrefs": [{"db_source": "MedGen", "db_id": "C1"}],
                    }
                ],
            },
        },
    }
}

GDC_PAGE = {
    "data": {
        "pagination": {"total": 2},
        "hits": [
            {
                "ssm": {
                    "ssm_id": "ssm-1",
                    "genomic_dna_change": "chr7:g.55191822T>G",
                    "consequence": [
                        {
                            "transcript": {
                                "transcript_id": "ENST00000275493",
                                "is_canonical": True,
                                "gene": {"symbol": "EGFR"},
                                "annotation": {"hgvsp_short": "p.L858R"},
                            }
                        }
                    ],
                },
                "case": {
                    "case_id": "case-1",
                    "submitter_id": "TCGA-01",
                    "project": {"project_id": "TCGA-LUAD"},
                },
            },
            {
                "ssm": {
                    "ssm_id": "ssm-2",
                    "consequence": [
                        {
                            "transcript": {
                                "is_canonical": True,
                                "gene": {"symbol": "EGFR"},
                                "annotation": {"hgvsp_short": "p.C240Y"},
                            }
                        }
                    ],
                },
                "case": {
                    "case_id": "case-2",
                    "submitter_id": "TCGA-02",
                    "project": {"project_id": "TCGA-GBM"},
                },
            },
        ],
    }
}


class NormalizationTests(unittest.TestCase):
    def test_normalizes_one_and_three_letter_hgvs(self):
        self.assertEqual(normalize_protein_change("p.Leu858Arg"), "L858R")
        self.assertEqual(normalize_protein_change("p.L858R"), "L858R")
        self.assertEqual(normalize_protein_change("p.Glu746_Ala750del"), "E746_A750del")

    def test_rejects_unknown_change(self):
        self.assertIsNone(normalize_protein_change("p.?"))


class EgfrPipelineTests(unittest.TestCase):
    def test_egfr_candidate_overlap(self):
        clinvar = NcbiClinVarClient(
            http=StubHttp(get_responses=[CLINVAR_SEARCH, CLINVAR_SUMMARY]),
            page_size=500,
        )
        gdc = GdcClient(
            http=StubHttp(get_responses=[GDC_STATUS], post_responses=[GDC_PAGE])
        )

        result = find_shared_mutations(
            "egfr", clinvar_client=clinvar, gdc_client=gdc
        )

        self.assertEqual(result.gene, "EGFR")
        self.assertEqual(result.rare_variant_count, 1)
        self.assertEqual(result.cancer_observation_count, 2)
        self.assertEqual(result.shared_mutation_count, 1)
        self.assertEqual(result.gdc_data_release, "Data Release fixture")
        self.assertEqual(len(result.mutation_plot_data), 2)
        plotted = {item.protein_change: item for item in result.mutation_plot_data}
        self.assertTrue(plotted["L858R"].in_rare_disease)
        self.assertTrue(plotted["L858R"].in_oncology)
        self.assertEqual(
            plotted["L858R"].project_patient_counts,
            (("TCGA-LUAD", 1),),
        )
        shared = result.shared_mutations[0]
        self.assertEqual(shared.protein_change, "L858R")
        self.assertEqual(shared.cancer_patient_count, 1)
        self.assertEqual(shared.tcga_projects, ("TCGA-LUAD",))
        self.assertEqual(shared.orphanet_codes, ("999999",))
        self.assertEqual(shared.genomic_match_count, 1)
        self.assertEqual(shared.match_level, "confirmed_genomic_and_protein")

    def test_minimum_patient_threshold(self):
        clinvar = NcbiClinVarClient(
            http=StubHttp(get_responses=[CLINVAR_SEARCH, CLINVAR_SUMMARY])
        )
        gdc = GdcClient(
            http=StubHttp(get_responses=[GDC_STATUS], post_responses=[GDC_PAGE])
        )
        result = find_shared_mutations(
            "EGFR",
            min_cancer_patients=2,
            clinvar_client=clinvar,
            gdc_client=gdc,
        )
        self.assertEqual(result.shared_mutation_count, 0)
        self.assertEqual(venn_counts(result).shared, 0)
        self.assertEqual(
            [
                (item.protein_change, item.in_rare_disease, item.in_oncology)
                for item in result.mutation_plot_data
            ],
            [("L858R", True, False)],
        )

    def test_can_exclude_a_generic_condition(self):
        clinvar = NcbiClinVarClient(
            http=StubHttp(get_responses=[CLINVAR_SEARCH, CLINVAR_SUMMARY])
        )
        gdc = GdcClient(
            http=StubHttp(get_responses=[GDC_STATUS], post_responses=[GDC_PAGE])
        )
        result = find_shared_mutations(
            "EGFR",
            exclude_conditions=["EGFR-associated rare disorder"],
            clinvar_client=clinvar,
            gdc_client=gdc,
        )
        self.assertEqual(result.rare_variant_count, 0)
        self.assertEqual(result.shared_mutation_count, 0)

    def test_partial_genomic_confirmation_is_not_propagated(self):
        rare = [
            RareVariant(
                gene="EGFR",
                protein_change="L858R",
                raw_protein_change="p.Leu858Arg",
                clinvar_id="101",
                clinvar_accession="VCV101",
                diseases=("Rare A",),
                orphanet_codes=("1",),
                condition_orphanet_pairs=(("Rare A", "1"),),
                germline_classification="Pathogenic",
                molecular_consequences=("missense variant",),
                canonical_spdi="NC_000007.14:55191821:T:G",
            ),
            RareVariant(
                gene="EGFR",
                protein_change="L858R",
                raw_protein_change="p.Leu858Arg",
                clinvar_id="202",
                clinvar_accession="VCV202",
                diseases=("Rare B",),
                orphanet_codes=("2",),
                condition_orphanet_pairs=(("Rare B", "2"),),
                germline_classification="Pathogenic",
                molecular_consequences=("missense variant",),
                canonical_spdi="NC_000007.14:999:T:G",
            ),
        ]
        cancer = [
            CancerObservation(
                gene="EGFR",
                protein_change="L858R",
                raw_protein_change="p.L858R",
                ssm_id="ssm-1",
                case_id="case-1",
                patient_id="TCGA-01",
                project_id="TCGA-LUAD",
                transcript_id="ENST00000275493",
                genomic_dna_change="chr7:g.55191822T>G",
            )
        ]

        result = _aggregate_matches(
            "EGFR", rare, cancer, min_cancer_patients=1
        )[0]
        self.assertEqual(result.genomic_match_count, 1)
        self.assertEqual(result.match_level, "mixed_genomic_and_candidate_protein")
        self.assertIn("VCV101|ssm-1", result.genomic_confirmations[0])
        self.assertNotIn("VCV202", result.genomic_confirmations[0])


class ClientRobustnessTests(unittest.TestCase):
    def test_clinvar_history_pagination(self):
        first = {
            "result": {
                "uids": ["101"],
                "101": CLINVAR_SUMMARY["result"]["101"],
            }
        }
        second = {
            "result": {
                "uids": ["102"],
                "102": CLINVAR_SUMMARY["result"]["102"],
            }
        }
        client = NcbiClinVarClient(
            http=StubHttp(get_responses=[CLINVAR_SEARCH, first, second]),
            page_size=1,
        )
        variants = client.fetch_rare_variants("EGFR")
        self.assertEqual([variant.clinvar_id for variant in variants], ["101"])

    def test_gdc_occurrence_pagination(self):
        first = {
            "data": {"pagination": {"total": 2}, "hits": [GDC_PAGE["data"]["hits"][0]]}
        }
        second = {
            "data": {"pagination": {"total": 2}, "hits": [GDC_PAGE["data"]["hits"][1]]}
        }
        client = GdcClient(
            http=StubHttp(
                get_responses=[GDC_STATUS],
                post_responses=[first, second],
            ),
            page_size=1,
        )
        observations = client.fetch_cancer_observations("EGFR")
        self.assertEqual(len(observations), 2)

    def test_api_key_is_redacted_from_http_errors(self):
        client = JsonHttpClient(user_agent="test", retries=0)
        url = "https://example.test/data?api_key=super-secret&gene=EGFR"
        error = urllib.error.HTTPError(url, 500, "failure", None, None)
        with mock.patch("urllib.request.urlopen", side_effect=error):
            with self.assertRaises(ApiError) as raised:
                client._request(url)
        self.assertNotIn("super-secret", str(raised.exception))
        self.assertIn("REDACTED", str(raised.exception))


if __name__ == "__main__":
    unittest.main()
