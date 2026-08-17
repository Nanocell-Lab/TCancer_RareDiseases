from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from rare_oncology.clients import NcbiClinVarClient
from rare_oncology.analysis import find_shared_mutations
from rare_oncology.evidence import (
    ClinVarEvidenceProvider,
    EvidenceContext,
    GnomadEvidenceProvider,
    LocalTableEvidenceProvider,
    OrphadataXmlEvidenceProvider,
    PanelAppEvidenceProvider,
    collect_evidence,
    spdi_identity,
)
from rare_oncology.export import write_evidence_csv, write_json
from rare_oncology.models import (
    AnalysisResult,
    ClinVarConditionAssertion,
    EvidenceRecord,
    RareVariant,
)


class StubHttp:
    def __init__(self, *, get_responses=(), post_responses=(), text_responses=()):
        self.get_responses = iter(get_responses)
        self.post_responses = iter(post_responses)
        self.text_responses = iter(text_responses)
        self.posts = []

    def get_json(self, _url, _params):
        return next(self.get_responses)

    def post_json(self, url, payload):
        self.posts.append((url, payload))
        return next(self.post_responses)

    def get_text(self, _url, _params):
        return next(self.text_responses)


def rare_variant() -> RareVariant:
    return RareVariant(
        gene="EGFR",
        protein_change="L858R",
        raw_protein_change="p.Leu858Arg",
        clinvar_id="101",
        clinvar_accession="VCV000000101.2",
        diseases=("Rare EGFR disorder",),
        orphanet_codes=("999",),
        condition_orphanet_pairs=(("Rare EGFR disorder", "999"),),
        germline_classification="Pathogenic",
        molecular_consequences=("missense variant",),
        canonical_spdi="NC_000007.14:55191821:T:G",
    )


class IdentityTests(unittest.TestCase):
    def test_spdi_becomes_exact_grch38_identity(self):
        identity = spdi_identity("NC_000007.14:55191821:T:G")
        self.assertEqual(identity["gnomad_variant_id"], "7-55191822-T-G")
        self.assertEqual(identity["variant_key"], "GRCh38:7:55191822:T>G")

    def test_unknown_or_non_literal_spdi_is_not_guessed(self):
        self.assertEqual(spdi_identity("NC_000007.14:1:3:G"), {})
        self.assertEqual(spdi_identity("NW_123.1:1:A:G"), {})

    def test_indel_or_ambiguous_spdi_is_not_sent_to_gnomad(self):
        deletion = spdi_identity("NC_000007.14:55191821:TC:T")
        self.assertEqual(deletion["variant_key"], "GRCh38:7:55191822:TC>T")
        self.assertNotIn("gnomad_variant_id", deletion)
        self.assertEqual(spdi_identity("NC_000007.14:1:N:G"), {})


class ClinVarDetailTests(unittest.TestCase):
    def test_vcv_xml_enriches_rcv_scv_and_review_provenance(self):
        search = {
            "esearchresult": {"count": "1", "querykey": "1", "webenv": "x"}
        }
        summary = {
            "result": {
                "uids": ["101"],
                "101": {
                    "uid": "101",
                    "accession_version": "VCV000000101.2",
                    "genes": [{"symbol": "EGFR"}],
                    "protein_change": "L858R",
                    "variation_set": [
                        {"canonical_spdi": "NC_000007.14:55191821:T:G"}
                    ],
                    "germline_classification": {
                        "description": "Pathogenic",
                        "trait_set": [
                            {
                                "trait_name": "Rare EGFR disorder",
                                "trait_xrefs": [
                                    {"db_source": "Orphanet", "db_id": "999"}
                                ],
                            }
                        ],
                    },
                },
            }
        }
        xml = """<ClinVarResult-Set>
          <VariationArchive VariationID="101" VariationType="single nucleotide variant"
             NumberOfSubmissions="2" NumberOfSubmitters="2" DateLastUpdated="2026-01-02">
            <ClassifiedRecord>
              <SimpleAllele><HGVSlist><HGVS Type="coding"><NucleotideExpression
                sequenceAccessionVersion="NM_005228.5" change="c.2573T&gt;G" />
              <ProteinExpression sequenceAccessionVersion="NP_005219.2"
                change="p.Leu858Arg" /></HGVS></HGVSlist></SimpleAllele>
              <RCVList>
                <RCVAccession Accession="RCV0001" Version="3">
                  <ClassifiedConditionList TraitSetID="1" />
                  <RCVClassifications><GermlineClassification>
                    <ReviewStatus>reviewed by expert panel</ReviewStatus>
                    <Description DateLastEvaluated="2025-12-01">Pathogenic</Description>
                  </GermlineClassification></RCVClassifications>
                </RCVAccession>
                <RCVAccession Accession="RCV0002" Version="1">
                  <ClassifiedConditionList TraitSetID="2" />
                  <RCVClassifications><GermlineClassification>
                    <ReviewStatus>criteria provided, single submitter</ReviewStatus>
                    <Description DateLastEvaluated="2024-01-01">Uncertain significance</Description>
                  </GermlineClassification></RCVClassifications>
                </RCVAccession>
              </RCVList>
              <Classifications>
                <ConditionList>
                  <TraitSet ID="1"><Trait><Name><ElementValue Type="Preferred">Disease A</ElementValue></Name>
                    <XRef DB="Orphanet" ID="1" /></Trait></TraitSet>
                  <TraitSet ID="2"><Trait><Name><ElementValue Type="Preferred">Disease B</ElementValue></Name>
                    <XRef DB="Orphanet" ID="2" /></Trait></TraitSet>
                </ConditionList>
                <GermlineClassification>
                <ReviewStatus>reviewed by expert panel</ReviewStatus>
                <Description DateLastEvaluated="2025-12-01">Pathogenic</Description>
              </GermlineClassification></Classifications>
              <ClinicalAssertionList><ClinicalAssertion><ClinVarAccession Type="SCV"
                Accession="SCV0001" Version="4" SubmitterName="Expert Lab" />
              </ClinicalAssertion></ClinicalAssertionList>
            </ClassifiedRecord>
          </VariationArchive>
        </ClinVarResult-Set>"""
        client = NcbiClinVarClient(
            http=StubHttp(
                get_responses=[search, summary], text_responses=[xml]
            )
        )
        variant = client.fetch_rare_variants("EGFR")[0]
        self.assertEqual(variant.rcv_accessions, ("RCV0001.3", "RCV0002.1"))
        self.assertEqual(variant.scv_accessions, ("SCV0001.4",))
        self.assertEqual(variant.submitters, ("Expert Lab",))
        self.assertEqual(variant.review_status, "reviewed by expert panel")
        self.assertEqual(variant.number_of_submissions, 2)
        self.assertIn("NM_005228.5:c.2573T>G", variant.hgvs_expressions)
        self.assertEqual(
            [(item.disease, item.classification) for item in variant.condition_assertions],
            [("Disease A", "Pathogenic"), ("Disease B", "Uncertain significance")],
        )

    def test_pipeline_filters_condition_specific_rcv_classification(self):
        variant = rare_variant()
        variant = RareVariant(
            **{
                **variant.__dict__,
                "diseases": ("Disease A", "Disease B"),
                "orphanet_codes": ("1", "2"),
                "condition_orphanet_pairs": (("Disease A", "1"), ("Disease B", "2")),
                "condition_assertions": (
                    ClinVarConditionAssertion(
                        "Disease A", ("1",), "RCV1.1", "Pathogenic", "expert panel", None
                    ),
                    ClinVarConditionAssertion(
                        "Disease B", ("2",), "RCV2.1", "Uncertain significance", "single submitter", None
                    ),
                ),
            }
        )

        class ClinVarStub:
            include_details = True

            def fetch_rare_variants(self, _gene, *, max_records=None):
                return [variant]

        class GdcStub:
            release_metadata = {}

            def fetch_cancer_observations(self, _gene, *, program, projects):
                return []

        result = find_shared_mutations(
            "EGFR", clinvar_client=ClinVarStub(), gdc_client=GdcStub()
        )
        disease_records = [
            record
            for record in result.evidence_records
            if record.source == "clinvar" and record.disease_id
        ]
        self.assertEqual(result.rare_variant_count, 1)
        self.assertEqual(
            [(record.disease_id, record.classification) for record in disease_records],
            [("ORPHA:1", "Pathogenic")],
        )

    def test_queried_details_never_fall_back_to_aggregate_classification(self):
        variant = rare_variant()

        class ClinVarStub:
            include_details = True

            def fetch_rare_variants(self, _gene, *, max_records=None):
                return [variant]

        class GdcStub:
            release_metadata = {}

            def fetch_cancer_observations(self, _gene, *, program, projects):
                return []

        result = find_shared_mutations(
            "EGFR", clinvar_client=ClinVarStub(), gdc_client=GdcStub()
        )
        status = next(
            item for item in result.source_statuses if item.source == "clinvar"
        )
        self.assertEqual(result.rare_variant_count, 0)
        self.assertEqual(status.status, "partial")
        self.assertIn("lack a resolvable RCV-Orphanet relationship", status.message)


class PublicProviderTests(unittest.TestCase):
    def setUp(self):
        self.context = EvidenceContext("EGFR", (rare_variant(),), ())

    def test_gnomad_maps_joint_counts_to_frequency(self):
        http = StubHttp(
            post_responses=[
                {
                    "data": {
                        "v0": {
                            "variant_id": "7-55191822-T-G",
                            "reference_genome": "GRCh38",
                            "chrom": "7",
                            "pos": 55191822,
                            "ref": "T",
                            "alt": "G",
                            "joint": {"ac": 2, "an": 100000},
                            "exome": {"ac": 2, "an": 90000},
                            "genome": {"ac": 0, "an": 10000},
                        }
                    }
                }
            ]
        )
        records = GnomadEvidenceProvider(http=http).collect(self.context)
        self.assertEqual(len(records), 1)
        self.assertAlmostEqual(records[0].population_frequency, 0.00002)
        self.assertIn("7-55191822-T-G", http.posts[0][1]["query"])

    def test_panelapp_keeps_gene_disease_evidence_separate(self):
        payload = {
            "results": [
                {
                    "gene_data": {"gene_symbol": "EGFR"},
                    "confidence_level": "3",
                    "phenotypes": ["Rare EGFR disorder"],
                    "mode_of_inheritance": "MONOALLELIC",
                    "evidence": ["Expert Review Green"],
                    "panel": {"id": 7, "name": "Rare panel", "version": "2.1"},
                }
            ]
        }
        records = PanelAppEvidenceProvider(
            http=StubHttp(get_responses=[payload])
        ).collect(self.context)
        self.assertEqual(records[0].evidence_layer, "gene_disease_validity")
        self.assertEqual(records[0].classification, "Green")
        self.assertIsNone(records[0].variant_key)

    def test_gnomad_not_found_is_an_explicit_absence_not_a_source_failure(self):
        provider = GnomadEvidenceProvider(
            http=StubHttp(
                post_responses=[
                    {
                        "data": {"v0": None},
                        "errors": [{"message": "Variant not found", "path": ["v0"]}],
                    }
                ]
            )
        )
        records, statuses = collect_evidence(self.context, (provider,))
        self.assertEqual(len(records), 1)
        self.assertIsNone(records[0].population_frequency)
        self.assertIn(("observed", "false"), records[0].details)
        self.assertEqual(statuses[0].status, "ok")
        self.assertIn("were not observed", statuses[0].message)

    def test_clinvar_submission_provenance_is_not_attributed_to_each_disease(self):
        variant = rare_variant()
        variant = RareVariant(
            **{
                **variant.__dict__,
                "condition_orphanet_pairs": (("Disease A", "1"), ("Disease B", "2")),
                "diseases": ("Disease A", "Disease B"),
                "orphanet_codes": ("1", "2"),
                "rcv_accessions": ("RCV1.1", "RCV2.1"),
                "scv_accessions": ("SCV1.1", "SCV2.1"),
                "submitters": ("Lab A", "Lab B"),
            }
        )
        records = ClinVarEvidenceProvider(details_queried=True).collect(
            EvidenceContext("EGFR", (variant,), ())
        )
        disease_records = [record for record in records if record.disease_id]
        provenance = [record for record in records if not record.disease_id]
        self.assertEqual(len(disease_records), 2)
        self.assertEqual(len(provenance), 1)
        self.assertNotIn("rcv_accessions", dict(disease_records[0].details))
        self.assertEqual(
            dict(provenance[0].details)["scope"],
            "variant_level_not_condition_attributed",
        )

    def test_gnomad_provider_resets_state_and_accumulates_warnings(self):
        base = rare_variant()
        deletion = RareVariant(
            **{**base.__dict__, "clinvar_id": "D", "canonical_spdi": "NC_000007.14:10:TC:T"}
        )
        second_snv = RareVariant(
            **{**base.__dict__, "clinvar_id": "S2", "canonical_spdi": "NC_000007.14:20:A:C"}
        )
        provider = GnomadEvidenceProvider(
            http=StubHttp(
                post_responses=[
                    {
                        "data": {"v0": None},
                        "errors": [{"message": "Variant not found", "path": ["v0"]}],
                    }
                ]
            ),
            max_variants=1,
        )
        provider.collect(EvidenceContext("EGFR", (base, second_snv, deletion), ()))
        self.assertIn("1 alleles skipped", provider.last_message)
        self.assertIn("Queried 1 of 2", provider.last_message)
        self.assertIn("1 of 1 alleles were not observed", provider.last_message)
        provider.collect(EvidenceContext("EGFR", (), ()))
        self.assertIsNone(provider.last_message)


class FileProviderTests(unittest.TestCase):
    def test_local_table_filters_gene_and_deduplicates_exact_records(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "erepo.tsv"
            path.write_text(
                "gene\trecord_id\thgvsp\tdisease_id\tclassification\n"
                "EGFR\tE1\tp.L858R\tMONDO:1\tPathogenic\n"
                "EGFR\tE1\tp.L858R\tMONDO:1\tPathogenic\n"
                "KRAS\tK1\tp.G12D\tMONDO:2\tPathogenic\n",
                encoding="utf-8",
            )
            provider = LocalTableEvidenceProvider("clingen", path)
            records, statuses = collect_evidence(
                EvidenceContext("EGFR", (), ()), (provider,)
            )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].protein_change, "L858R")
        self.assertEqual(statuses[0].record_count, 2)

    def test_orphadata_xml_returns_gene_disease_association(self):
        xml = """<JDBOR date="2026-07-01" version="test">
          <DisorderList><Disorder><OrphaCode>337</OrphaCode>
          <Name>Fibrodysplasia ossificans progressiva</Name>
          <DisorderGeneAssociationList><DisorderGeneAssociation>
          <SourceOfValidation>123[PMID]</SourceOfValidation>
          <Gene><Symbol>ACVR1</Symbol><ExternalReferenceList>
          <ExternalReference><Source>HGNC</Source><Reference>171</Reference></ExternalReference>
          </ExternalReferenceList></Gene>
          <DisorderGeneAssociationType><Name>Disease-causing germline mutation(s) in</Name></DisorderGeneAssociationType>
          <DisorderGeneAssociationStatus><Name>Assessed</Name></DisorderGeneAssociationStatus>
          </DisorderGeneAssociation></DisorderGeneAssociationList>
          </Disorder></DisorderList></JDBOR>"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "orpha.xml"
            path.write_text(xml, encoding="utf-8")
            provider = OrphadataXmlEvidenceProvider(path)
            records = provider.collect(EvidenceContext("ACVR1", (), ()))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].disease_id, "ORPHA:337")
        self.assertEqual(records[0].evidence_level, "Assessed")
        self.assertEqual(records[0].license, "CC BY 4.0")

    def test_optional_failure_is_visible_and_strict_mode_raises(self):
        provider = LocalTableEvidenceProvider("lovd", "/missing/evidence.tsv")
        records, statuses = collect_evidence(
            EvidenceContext("EGFR", (), ()), (provider,)
        )
        self.assertEqual(records, ())
        self.assertEqual(statuses[0].status, "error")
        with self.assertRaises(ValueError):
            collect_evidence(
                EvidenceContext("EGFR", (), ()), (provider,), strict=True
            )

    def test_local_table_requires_an_explicit_gene(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "unsafe.csv"
            path.write_text(
                "record_id,disease\nKRAS-1,KRAS syndrome\n", encoding="utf-8"
            )
            provider = LocalTableEvidenceProvider("lovd", path)
            records, statuses = collect_evidence(
                EvidenceContext("ACVR1", (), ()), (provider,)
            )
        self.assertEqual(records, ())
        self.assertEqual(statuses[0].status, "error")
        self.assertIn("has no gene", statuses[0].message)

    def test_invalid_population_frequency_is_rejected(self):
        for value in ("nan", "inf", "-0.1", "1.1"):
            with self.subTest(value=value), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "frequency.csv"
                path.write_text(
                    f"gene,record_id,af\nEGFR,E1,{value}\n", encoding="utf-8"
                )
                records, statuses = collect_evidence(
                    EvidenceContext("EGFR", (), ()),
                    (LocalTableEvidenceProvider("gnomad_export", path),),
                )
                self.assertEqual(records, ())
                self.assertEqual(statuses[0].status, "error")


class ExportSafetyTests(unittest.TestCase):
    def test_csv_formula_is_escaped_and_json_rejects_nan(self):
        record = EvidenceRecord(
            source="local",
            source_record_id="=HYPERLINK(\"https://example.invalid\")",
            evidence_layer="variant_assertion",
            gene="EGFR",
            retrieved_at="2026-01-01T00:00:00+00:00",
            population_frequency=float("nan"),
        )
        result = AnalysisResult(
            gene="EGFR",
            generated_at="2026-01-01T00:00:00+00:00",
            gdc_program="TCGA",
            gdc_projects=(),
            min_cancer_patients=1,
            excluded_conditions=(),
            clinvar_record_limit=None,
            allowed_germline_classifications=None,
            gdc_data_release=None,
            rare_variant_count=0,
            cancer_observation_count=0,
            shared_mutation_count=0,
            shared_mutations=(),
            mutation_plot_data=(),
            evidence_record_count=1,
            evidence_records=(record,),
        )
        with tempfile.TemporaryDirectory() as directory:
            csv_path = write_evidence_csv(result, Path(directory) / "evidence.csv")
            self.assertIn("'=HYPERLINK", csv_path.read_text(encoding="utf-8"))
            with self.assertRaises(ValueError):
                write_json(result, Path(directory) / "report.json")


if __name__ == "__main__":
    unittest.main()
