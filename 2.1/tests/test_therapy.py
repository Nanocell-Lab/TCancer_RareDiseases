from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from rare_oncology.chemical_metadata import ChemicalMetadata
from rare_oncology.cli import _build_therapy_providers, build_parser
from rare_oncology.export import write_therapy_csv, write_therapy_json
from rare_oncology.models import AnalysisResult, SharedMutation, TherapyAssociation
from rare_oncology.therapy import (
    REGULATORY_THERAPY,
    TRIAL_AND_TARGET_DISCOVERY,
    VARIANT_THERAPY_EVIDENCE,
    CivicTherapyProvider,
    ClinicalTrialsTherapyProvider,
    OncoKBTherapyProvider,
    OpenTargetsTherapyProvider,
    RegulatoryTableTherapyProvider,
    TherapyContext,
    collect_therapy_associations,
    find_therapy_associations,
)


class StubHttp:
    def __init__(self, *, gets=(), posts=()):
        self.gets = iter(gets)
        self.posts = iter(posts)
        self.get_calls = []
        self.get_headers = []
        self.post_calls = []

    def get_json(self, url, params, *, headers=None):
        self.get_calls.append((url, params))
        self.get_headers.append(headers)
        return next(self.gets)

    def post_json(self, url, payload):
        self.post_calls.append((url, payload))
        return next(self.posts)


class StubChemicalResolver:
    def __init__(self, metadata=(), error=None):
        self.metadata = {item.chembl_id: item for item in metadata}
        self.error = error
        self.calls = []

    def cache_identity(self):
        return {"source": "chemical-fixture"}

    def resolve(self, chembl_ids):
        self.calls.append(chembl_ids)
        if self.error is not None:
            raise self.error
        return {
            chembl_id: self.metadata[chembl_id]
            for chembl_id in chembl_ids
            if chembl_id in self.metadata
        }


def analysis_result() -> AnalysisResult:
    shared = SharedMutation(
        gene="EGFR",
        protein_change="L858R",
        diseases=("Rare EGFR disorder",),
        orphanet_codes=("999",),
        clinvar_ids=("101",),
        clinvar_accessions=("VCV000000101.2",),
        germline_classifications=("Pathogenic",),
        canonical_spdis=("NC_000007.14:55191821:T:G",),
        tcga_projects=("TCGA-LUAD",),
        cancer_patient_count=2,
        cancer_event_count=2,
        genomic_match_count=1,
        genomic_confirmations=("GRCh38:7:55191822:T>G",),
        ssm_ids=("ssm-1",),
        match_level="confirmed",
    )
    return AnalysisResult(
        gene="EGFR",
        generated_at="2026-08-11T00:00:00+00:00",
        gdc_program="TCGA",
        gdc_projects=("TCGA-LUAD",),
        min_cancer_patients=1,
        excluded_conditions=(),
        clinvar_record_limit=None,
        allowed_germline_classifications=("Pathogenic",),
        gdc_data_release="42",
        rare_variant_count=1,
        cancer_observation_count=2,
        shared_mutation_count=1,
        shared_mutations=(shared,),
        mutation_plot_data=(),
    )


class ContextTests(unittest.TestCase):
    def test_context_uses_only_shared_mutations(self):
        context = TherapyContext.from_analysis(analysis_result())
        self.assertEqual(context.gene, "EGFR")
        self.assertEqual(context.protein_changes, ("L858R",))
        self.assertEqual(context.disease_ids, ("ORPHA:999",))


class CivicTests(unittest.TestCase):
    def test_accepted_predictive_evidence_is_variant_specific(self):
        variants = {
            "data": {
                "gene": {
                    "variants": {
                        "pageInfo": {"hasNextPage": False, "endCursor": "x"},
                        "nodes": [
                            {"id": 1, "name": "L858R"},
                            {"id": 2, "name": "T790M"},
                        ],
                    }
                }
            }
        }
        evidence = {
            "data": {
                "evidenceItems": {
                    "pageInfo": {"hasNextPage": False, "endCursor": None},
                    "nodes": [
                        {
                            "id": 900,
                            "status": "ACCEPTED",
                            "evidenceType": "PREDICTIVE",
                            "evidenceLevel": "A",
                            "evidenceDirection": "SUPPORTS",
                            "significance": "SENSITIVITYRESPONSE",
                            "evidenceRating": 5,
                            "therapyInteractionType": None,
                            "description": "Accepted evidence",
                            "link": "/evidence/900",
                            "molecularProfile": {
                                "id": 20,
                                "name": "EGFR L858R",
                                "isMultiVariant": False,
                            },
                            "disease": {
                                "id": 10,
                                "name": "Lung Adenocarcinoma",
                                "doid": "3910",
                            },
                            "therapies": [
                                {
                                    "id": 30,
                                    "name": "Exampletinib",
                                    "ncitId": "C123",
                                    "myChemInfo": {
                                        "chemblId": "CHEMBL123",
                                        "chemblMoleculeType": "Small molecule",
                                    },
                                }
                            ],
                            "source": {
                                "citationId": "123456",
                                "sourceType": "PUBMED",
                            },
                        }
                    ],
                }
            }
        }
        http = StubHttp(posts=[variants, evidence])
        resolver = StubChemicalResolver(
            [ChemicalMetadata("CHEMBL123", "Small molecule", "CCO")]
        )
        records = CivicTherapyProvider(
            http=http,
            chemical_resolver=resolver,
        ).collect(
            TherapyContext.from_analysis(analysis_result())
        )
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record.layer, VARIANT_THERAPY_EVIDENCE)
        self.assertEqual(record.association_scope, "candidate_protein_variant")
        self.assertEqual(record.protein_change, "L858R")
        self.assertEqual(record.therapy_name, "Exampletinib")
        self.assertEqual(record.effect, "sensitivity_response")
        self.assertEqual(record.therapy_status, "curated_predictive_evidence")
        self.assertEqual(record.therapy_type, "small_molecule")
        self.assertEqual(record.chembl_ids, ("CHEMBL123",))
        self.assertEqual(record.canonical_smiles, "CCO")
        self.assertEqual(resolver.calls, [("CHEMBL123",)])
        requested_ids = [
            call[1]["variables"].get("variantId") for call in http.post_calls
        ]
        self.assertEqual(requested_ids, [None, 1])

    def test_graphql_errors_are_reported_as_source_error(self):
        provider = CivicTherapyProvider(
            http=StubHttp(posts=[{"errors": [{"message": "schema changed"}]}])
        )
        records, statuses = collect_therapy_associations(
            TherapyContext.from_analysis(analysis_result()), [provider]
        )
        self.assertEqual(records, ())
        self.assertEqual(statuses[0].status, "error")
        self.assertIsNotNone(statuses[0].message)
        self.assertIn("schema changed", statuses[0].message or "")

    def test_missing_graphql_cursor_is_a_visible_error(self):
        provider = CivicTherapyProvider(
            http=StubHttp(
                posts=[
                    {
                        "data": {
                            "gene": {
                                "variants": {
                                    "pageInfo": {
                                        "hasNextPage": True,
                                        "endCursor": None,
                                    },
                                    "nodes": [],
                                }
                            }
                        }
                    }
                ]
            )
        )
        records, statuses = collect_therapy_associations(
            TherapyContext.from_analysis(analysis_result()), [provider]
        )
        self.assertEqual(records, ())
        self.assertEqual(statuses[0].status, "error")
        self.assertIn("missing cursor", statuses[0].message or "")


class OncoKBTests(unittest.TestCase):
    def test_variant_treatments_preserve_level_tumor_effect_and_bearer_auth(self):
        response = {
            "dataVersion": "v7.4",
            "lastUpdate": "07/31/2026",
            "oncogenic": "Oncogenic",
            "highestSensitiveLevel": "LEVEL_1",
            "highestResistanceLevel": "LEVEL_R1",
            "treatments": [
                {
                    "drugs": [
                        {"ncitCode": "C1", "drugName": "Alpha"},
                        {"ncitCode": "C2", "drugName": "Beta"},
                    ],
                    "level": "LEVEL_1",
                    "fdaLevel": "LEVEL_Fda2",
                    "levelAssociatedCancerType": {
                        "id": 78,
                        "code": "",
                        "name": "",
                        "mainType": {"name": "Non-Small Cell Lung Cancer"},
                    },
                    "pmids": ["123", "456"],
                    "description": "licensed text must not be redistributed",
                },
                {
                    "drugs": [{"ncitCode": "C3", "drugName": "Gamma"}],
                    "level": "LEVEL_R1",
                    "levelAssociatedCancerType": {
                        "code": "NSCLC",
                        "name": "Non-Small Cell Lung Cancer",
                    },
                    "pmids": ["789"],
                },
            ],
        }
        http = StubHttp(gets=[response])
        provider = OncoKBTherapyProvider("test-token", http=http)
        records = provider.collect(TherapyContext.from_analysis(analysis_result()))
        self.assertEqual(len(records), 2)
        by_level = {record.evidence_level: record for record in records}
        sensitive = by_level["LEVEL_1"]
        resistant = by_level["LEVEL_R1"]
        self.assertEqual(sensitive.layer, VARIANT_THERAPY_EVIDENCE)
        self.assertEqual(sensitive.association_scope, "oncokb_protein_variant")
        self.assertEqual(sensitive.therapy_components, ("Alpha", "Beta"))
        self.assertEqual(sensitive.therapy_type, "combination_drug")
        self.assertEqual(sensitive.therapy_id, "NCIT:C1 + NCIT:C2")
        self.assertEqual(sensitive.disease_ids, ("ONCOTREE_ID:78",))
        self.assertEqual(
            sensitive.disease_names, ("Non-Small Cell Lung Cancer",)
        )
        self.assertEqual(sensitive.effect, "sensitivity_response")
        self.assertEqual(resistant.effect, "resistance")
        self.assertEqual(resistant.therapy_type, "drug")
        self.assertIsNone(sensitive.description)
        self.assertEqual(provider.source_version, "v7.4")
        self.assertEqual(
            http.get_headers,
            [{"Authorization": "Bearer test-token"}],
        )
        self.assertNotIn("test-token", json.dumps(provider.cache_identity()))

    def test_variant_limit_is_reported_as_partial(self):
        provider = OncoKBTherapyProvider(
            "test-token",
            http=StubHttp(gets=[{"dataVersion": "v7.4", "treatments": []}]),
            max_variants=1,
        )
        records, statuses = collect_therapy_associations(
            TherapyContext("EGFR", ("L858R", "T790M"), (), ()), [provider]
        )
        self.assertEqual(records, ())
        self.assertEqual(statuses[0].status, "partial")
        self.assertIn("1/2", statuses[0].message or "")


class ClinicalTrialsTests(unittest.TestCase):
    def test_trial_with_variant_is_not_called_approved_or_effective(self):
        response = {
            "studies": [
                {
                    "protocolSection": {
                        "identificationModule": {
                            "nctId": "NCT00000001",
                            "briefTitle": "Study for EGFR L858R tumors",
                        },
                        "statusModule": {"overallStatus": "RECRUITING"},
                        "designModule": {"phases": ["PHASE2"]},
                        "conditionsModule": {
                            "conditions": ["Lung Adenocarcinoma"]
                        },
                        "armsInterventionsModule": {
                            "interventions": [
                                {"type": "DRUG", "name": "Trialtinib"},
                                {"type": "PROCEDURE", "name": "Biopsy"},
                            ]
                        },
                        "eligibilityModule": {
                            "eligibilityCriteria": "EGFR p.L858R required"
                        },
                    }
                }
            ]
        }
        http = StubHttp(gets=[response, {"studies": []}])
        records = ClinicalTrialsTherapyProvider(http=http).collect(
            TherapyContext.from_analysis(analysis_result())
        )
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record.layer, TRIAL_AND_TARGET_DISCOVERY)
        self.assertEqual(record.association_scope, "trial_variant_mention")
        self.assertEqual(record.therapy_status, "clinical_trial")
        self.assertEqual(record.trial_status, "RECRUITING")
        self.assertEqual(record.therapy_type, "drug")
        self.assertNotEqual(record.therapy_status, "regulatory_approved")
        self.assertIn("eligibility", dict(record.details)["match_locations"])
        self.assertEqual(
            [call[1]["query.term"] for call in http.get_calls],
            ["L858R", "EGFR gene"],
        )

    def test_short_gene_symbol_requires_uppercase_biomedical_token(self):
        response = {
            "studies": [
                {
                    "protocolSection": {
                        "identificationModule": {
                            "nctId": "NCT00000002",
                            "briefTitle": "Primary endpoint was met",
                        },
                        "armsInterventionsModule": {
                            "interventions": [
                                {"type": "DRUG", "name": "Unrelated drug"}
                            ]
                        },
                    }
                }
            ]
        }
        records = ClinicalTrialsTherapyProvider(
            http=StubHttp(gets=[response])
        ).collect(TherapyContext("MET", (), (), ()))
        self.assertEqual(records, [])

    def test_placebo_is_excluded_and_active_comparator_is_labelled(self):
        response = {
            "studies": [
                {
                    "protocolSection": {
                        "identificationModule": {
                            "nctId": "NCT00000003",
                            "briefTitle": "EGFR randomized trial",
                        },
                        "armsInterventionsModule": {
                            "interventions": [
                                {"type": "DRUG", "name": "Study drug"},
                                {"type": "DRUG", "name": "Placebo"},
                                {"type": "DRUG", "name": "Standard therapy"},
                            ],
                            "armGroups": [
                                {
                                    "type": "EXPERIMENTAL",
                                    "interventionNames": ["Drug: Study drug"],
                                },
                                {
                                    "type": "PLACEBO_COMPARATOR",
                                    "interventionNames": ["Drug: Placebo"],
                                },
                                {
                                    "type": "ACTIVE_COMPARATOR",
                                    "interventionNames": ["Drug: Standard therapy"],
                                },
                            ],
                        },
                    }
                }
            ]
        }
        records = ClinicalTrialsTherapyProvider(
            http=StubHttp(gets=[response, {"studies": []}])
        ).collect(TherapyContext.from_analysis(analysis_result()))
        self.assertEqual(
            [(item.therapy_name, item.trial_role) for item in records],
            [
                ("Standard therapy", "active_comparator"),
                ("Study drug", "experimental"),
            ],
        )

    def test_repeated_pagination_token_is_a_visible_error(self):
        page = {
            "studies": [
                {
                    "protocolSection": {
                        "identificationModule": {
                            "nctId": "NCT00000004",
                            "briefTitle": "EGFR study",
                        }
                    }
                }
            ],
            "nextPageToken": "same",
        }
        provider = ClinicalTrialsTherapyProvider(
            http=StubHttp(gets=[page, page]), max_studies=3
        )
        records, statuses = collect_therapy_associations(
            TherapyContext.from_analysis(analysis_result()), [provider]
        )
        self.assertEqual(records, ())
        self.assertEqual(statuses[0].status, "error")
        self.assertIn("repeated", statuses[0].message or "")

    def test_raw_study_cap_is_reported_as_partial(self):
        response = {
            "studies": [
                {
                    "protocolSection": {
                        "identificationModule": {
                            "nctId": "NCT00000005",
                            "briefTitle": "MET inhibitor study",
                        },
                        "armsInterventionsModule": {
                            "interventions": [
                                {"type": "DRUG", "name": "Metinib"}
                            ]
                        },
                    }
                }
            ],
            "nextPageToken": "more",
        }
        provider = ClinicalTrialsTherapyProvider(
            http=StubHttp(gets=[response]), max_studies=1
        )
        records, statuses = collect_therapy_associations(
            TherapyContext("MET", (), (), ()), [provider]
        )
        self.assertEqual(len(records), 1)
        self.assertEqual(statuses[0].status, "partial")
        self.assertIn("truncated", statuses[0].message or "")


class OpenTargetsTests(unittest.TestCase):
    def test_approval_stage_remains_gene_target_discovery(self):
        search = {
            "data": {
                "search": {
                    "hits": [
                        {
                            "id": "ENSG00000146648",
                            "name": "EGFR",
                            "entity": "target",
                        }
                    ]
                }
            }
        }
        drugs = {
            "data": {
                "target": {
                    "id": "ENSG00000146648",
                    "drugAndClinicalCandidates": {
                        "rows": [
                            {
                                "id": "row-1",
                                "maxClinicalStage": "APPROVAL",
                                "drug": {
                                    "id": "CHEMBL1",
                                    "name": "Targetinib",
                                    "drugType": "Small molecule",
                                    "maximumClinicalStage": "APPROVAL",
                                },
                                "diseases": [
                                    {
                                        "diseaseFromSource": "lung cancer",
                                        "disease": {
                                            "id": "MONDO_0005233",
                                            "name": "lung carcinoma",
                                        },
                                    }
                                ],
                            }
                        ]
                    },
                }
            }
        }
        resolver = StubChemicalResolver(
            [ChemicalMetadata("CHEMBL1", "Small molecule", "C1=CC=CC=C1")]
        )
        records = OpenTargetsTherapyProvider(
            http=StubHttp(posts=[search, drugs]),
            chemical_resolver=resolver,
        ).collect(TherapyContext.from_analysis(analysis_result()))
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record.association_scope, "gene_target")
        self.assertEqual(record.therapy_status, "mechanism_target")
        self.assertEqual(record.evidence_level, "APPROVAL")
        self.assertEqual(record.therapy_type, "small_molecule")
        self.assertEqual(record.chembl_ids, ("CHEMBL1",))
        self.assertEqual(record.canonical_smiles, "C1=CC=CC=C1")

    def test_smiles_failure_preserves_therapy_and_marks_source_partial(self):
        search = {
            "data": {
                "search": {
                    "hits": [
                        {
                            "id": "ENSG00000146648",
                            "name": "EGFR",
                            "entity": "target",
                        }
                    ]
                }
            }
        }
        drugs = {
            "data": {
                "target": {
                    "id": "ENSG00000146648",
                    "drugAndClinicalCandidates": {
                        "rows": [
                            {
                                "id": "row-1",
                                "drug": {
                                    "id": "CHEMBL1",
                                    "name": "Targetinib",
                                    "drugType": "Small molecule",
                                },
                                "diseases": [],
                            }
                        ]
                    },
                }
            }
        }
        provider = OpenTargetsTherapyProvider(
            http=StubHttp(posts=[search, drugs]),
            chemical_resolver=StubChemicalResolver(error=OSError("offline")),
        )

        records, statuses = collect_therapy_associations(
            TherapyContext.from_analysis(analysis_result()), [provider]
        )

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].therapy_type, "small_molecule")
        self.assertIsNone(records[0].canonical_smiles)
        self.assertEqual(statuses[0].status, "partial")
        self.assertIn("SMILES", statuses[0].message or "")


class RegulatoryTests(unittest.TestCase):
    def test_official_table_preserves_jurisdiction_and_status(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fda.csv"
            path.write_text(
                "gene,record_id,drug,variant,indication,regulatory_status,label_url,"
                "therapy_type,chembl_id,canonical_smiles\n"
                "EGFR,NDA-1,Labeltinib,p.L858R,Lung cancer,regulatory_approved,"
                "https://example.test/label,Small molecule,chembl42,CCN\n",
                encoding="utf-8",
            )
            records = RegulatoryTableTherapyProvider(
                "fda", path, jurisdiction="US"
            ).collect(TherapyContext.from_analysis(analysis_result()))
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record.layer, REGULATORY_THERAPY)
        self.assertEqual(record.therapy_status, "regulatory_approved")
        self.assertEqual(record.jurisdiction, "US")
        self.assertEqual(record.protein_change, "L858R")
        self.assertEqual(record.therapy_type, "small_molecule")
        self.assertEqual(record.chembl_ids, ("CHEMBL42",))
        self.assertEqual(record.canonical_smiles, "CCN")

    def test_gene_is_required_per_regulatory_row(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fda.csv"
            path.write_text("record_id,drug\n1,UnsafeGuess\n", encoding="utf-8")
            provider = RegulatoryTableTherapyProvider("fda", path)
            records, statuses = collect_therapy_associations(
                TherapyContext.from_analysis(analysis_result()), [provider]
            )
        self.assertEqual(records, ())
        self.assertEqual(statuses[0].status, "error")
        self.assertIsNotNone(statuses[0].message)
        self.assertIn("has no gene", statuses[0].message or "")

    def test_same_application_keeps_distinct_indications(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fda.csv"
            path.write_text(
                "gene,record_id,drug,indication,regulatory_status\n"
                "EGFR,NDA-1,Labeltinib,Lung cancer,regulatory_approved\n"
                "EGFR,NDA-1,Labeltinib,Glioma,label_biomarker\n",
                encoding="utf-8",
            )
            records = RegulatoryTableTherapyProvider("fda", path).collect(
                TherapyContext.from_analysis(analysis_result())
            )
        self.assertEqual(len(records), 2)

    def test_unstructured_variant_is_preserved_as_biomarker(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fda.csv"
            path.write_text(
                "gene,record_id,drug,variant,regulatory_status\n"
                "EGFR,NDA-2,Labeltinib,EGFR exon 19 deletion,label_biomarker\n",
                encoding="utf-8",
            )
            record = RegulatoryTableTherapyProvider("fda", path).collect(
                TherapyContext.from_analysis(analysis_result())
            )[0]
        self.assertIsNone(record.protein_change)
        self.assertEqual(record.association_scope, "biomarker")
        self.assertEqual(record.molecular_profile, "EGFR exon 19 deletion")

    def test_missing_regulatory_status_is_explicitly_unknown(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fda.csv"
            path.write_text(
                "gene,record_id,drug\nEGFR,NDA-3,UnknownStatusDrug\n",
                encoding="utf-8",
            )
            record = RegulatoryTableTherapyProvider("fda", path).collect(
                TherapyContext.from_analysis(analysis_result())
            )[0]
        self.assertEqual(record.therapy_status, "regulatory_status_not_provided")


class PartitionAndExportTests(unittest.TestCase):
    def test_cli_exposes_opt_in_therapy_sources_and_regulatory_files(self):
        args = build_parser().parse_args(
            [
                "EGFR",
                "--therapy-source",
                "civic",
                "--regulatory-file",
                "fda=fda.csv",
            ]
        )
        self.assertEqual(args.therapy_source, ["civic"])
        self.assertEqual(args.regulatory_file, [("fda", Path("fda.csv"))])

    def test_cli_builds_oncokb_only_with_explicit_token(self):
        args = build_parser().parse_args(
            [
                "EGFR",
                "--therapy-source",
                "oncokb",
                "--oncokb-token",
                "test-token",
            ]
        )
        providers = _build_therapy_providers(args)
        self.assertEqual([provider.name for provider in providers], ["oncokb"])
        missing = build_parser().parse_args(
            ["EGFR", "--therapy-source", "oncokb", "--oncokb-token", ""]
        )
        with self.assertRaisesRegex(SystemExit, "ONCOKB_TOKEN"):
            _build_therapy_providers(missing)

    def test_regulatory_file_alone_does_not_enable_external_queries(self):
        args = build_parser().parse_args(
            ["EGFR", "--regulatory-file", "fda=fda.csv"]
        )
        providers = _build_therapy_providers(args)
        self.assertEqual([provider.name for provider in providers], ["fda"])

    def test_result_partitions_records_by_layer(self):
        class Provider:
            name = "fixture"
            source_version = "1"
            url = None

            def __init__(self, layer):
                self.layer = layer

            def collect(self, context):
                return [
                    TherapyAssociation(
                        layer=self.layer,
                        source=self.name,
                        source_record_id=self.layer,
                        gene=context.gene,
                        association_scope="gene_target",
                        therapy_name=self.layer,
                        therapy_components=(self.layer,),
                        therapy_status="fixture",
                        retrieved_at="2026-08-11T00:00:00+00:00",
                    )
                ]

        result = find_therapy_associations(
            analysis_result(),
            providers=[
                Provider(VARIANT_THERAPY_EVIDENCE),
                Provider(REGULATORY_THERAPY),
                Provider(TRIAL_AND_TARGET_DISCOVERY),
            ],
        )
        self.assertEqual(result.association_count, 3)
        self.assertEqual(len(result.variant_therapy_evidence), 1)
        self.assertEqual(len(result.regulatory_therapy), 1)
        self.assertEqual(len(result.trial_and_target_discovery), 1)

    def test_unknown_provider_layer_is_rejected_visibly(self):
        class Provider:
            name = "typo"
            layer = "regluatory_therapy"
            source_version = None
            url = None

            def collect(self, context):
                raise AssertionError("collect no debe ejecutarse")

        records, statuses = collect_therapy_associations(
            TherapyContext.from_analysis(analysis_result()), [Provider()]
        )
        self.assertEqual(records, ())
        self.assertEqual(statuses[0].status, "error")
        self.assertIn("unknown", statuses[0].message or "")

    def test_json_and_csv_exports_are_nonempty_and_safe(self):
        class Provider:
            name = "fixture"
            layer = REGULATORY_THERAPY
            source_version = "1"
            url = None

            def collect(self, context):
                return [
                    TherapyAssociation(
                        layer=self.layer,
                        source=self.name,
                        source_record_id="1",
                        gene=context.gene,
                        association_scope="biomarker",
                        therapy_name="\t=FORMULA()",
                        therapy_components=("\t=FORMULA()",),
                        therapy_status="label_biomarker",
                        retrieved_at="2026-08-11T00:00:00+00:00",
                        therapy_type="small_molecule",
                        chembl_ids=("CHEMBL42",),
                        canonical_smiles="CCN",
                    )
                ]

        result = find_therapy_associations(
            analysis_result(), providers=[Provider()]
        )
        with tempfile.TemporaryDirectory() as directory:
            json_path = write_therapy_json(result, Path(directory) / "therapy.json")
            csv_path = write_therapy_csv(result, Path(directory) / "therapy.csv")
            payload = json.loads(json_path.read_text(encoding="utf-8"))
            with csv_path.open(encoding="utf-8") as handle:
                row = next(csv.DictReader(handle))
        self.assertEqual(payload["association_count"], 1)
        self.assertEqual(
            payload["regulatory_therapy"][0]["canonical_smiles"], "CCN"
        )
        self.assertEqual(row["therapy_type"], "small_molecule")
        self.assertEqual(row["chembl_ids"], "CHEMBL42")
        self.assertEqual(row["canonical_smiles"], "CCN")
        self.assertTrue(row["therapy_name"].startswith("'\t="))


if __name__ == "__main__":
    unittest.main()
