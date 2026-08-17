from __future__ import annotations

import csv
import sqlite3
import tempfile
import unittest
from pathlib import Path

from rare_oncology.alphamissense import AlphaMissenseEvidenceProvider
from rare_oncology.evidence import EvidenceContext
from rare_oncology.export import write_csv, write_evidence_csv
from rare_oncology.models import (
    AnalysisResult,
    EvidenceRecord,
    MutationPlotRecord,
    RareVariant,
    SharedMutation,
)


class StubHttp:
    def __init__(self, responses=()):
        self.responses = iter(responses)
        self.posts: list[tuple[str, dict]] = []

    def post_json(self, url, payload):
        self.posts.append((url, payload))
        return next(self.responses)


def rare_variant(
    *,
    change: str = "L858R",
    spdi: str = "NC_000007.14:55191821:T:G",
) -> RareVariant:
    return RareVariant(
        gene="EGFR",
        protein_change=change,
        raw_protein_change=f"p.{change}",
        clinvar_id="101",
        clinvar_accession="VCV000000101.2",
        diseases=("Rare EGFR disorder",),
        orphanet_codes=("999",),
        condition_orphanet_pairs=(("Rare EGFR disorder", "999"),),
        germline_classification="Pathogenic",
        molecular_consequences=("missense variant",),
        canonical_spdi=spdi,
    )


def context(*variants: RareVariant) -> EvidenceContext:
    return EvidenceContext(
        gene="EGFR", rare_variants=tuple(variants), cancer_observations=()
    )


def vep_response(
    *,
    identifier: str = "am0",
    amino_acids: str = "L/R",
    score: float = 0.9968,
    prediction_class: str = "likely_pathogenic",
):
    return [
        {
            "id": identifier,
            "transcript_consequences": [
                {
                    "gene_symbol": "EGFR",
                    "transcript_id": "ENST00000275493",
                    "consequence_terms": ["missense_variant"],
                    "protein_start": 858,
                    "amino_acids": amino_acids,
                    "canonical": 1,
                    "mane_select": "NM_005228.5",
                    "alphamissense": {
                        "am_pathogenicity": score,
                        "am_class": prediction_class,
                    },
                },
                {
                    "gene_symbol": "OTHER",
                    "transcript_id": "ENST_OTHER",
                    "consequence_terms": ["missense_variant"],
                    "protein_start": 1,
                    "amino_acids": "A/V",
                    "alphamissense": {
                        "am_pathogenicity": 0.8,
                        "am_class": "likely_pathogenic",
                    },
                },
            ],
        }
    ]


class AlphaMissenseProviderTests(unittest.TestCase):
    def test_exact_genomic_and_protein_match_returns_score(self):
        http = StubHttp((vep_response(),))
        provider = AlphaMissenseEvidenceProvider(http=http)

        records = provider.collect(context(rare_variant()))

        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record.variant_key, "GRCh38:7:55191822:T>G")
        self.assertEqual(record.protein_change, "L858R")
        self.assertEqual(record.transcript, "ENST00000275493")
        self.assertEqual(record.prediction_score, 0.9968)
        self.assertEqual(record.prediction_class, "likely_pathogenic")
        self.assertEqual(record.evidence_level, "computational_prediction")
        self.assertIn("AlphaMissense=1", http.posts[0][0])
        self.assertEqual(
            http.posts[0][1]["variants"],
            ["7 55191822 am0 T G . . ."],
        )

    def test_sqlite_reuses_found_and_unavailable_variants(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "runs.sqlite3"
            first_http = StubHttp((vep_response(),))
            first = AlphaMissenseEvidenceProvider(
                database=database, http=first_http
            )
            first_records = first.collect(context(rare_variant()))

            cached_http = StubHttp()
            cached = AlphaMissenseEvidenceProvider(
                database=database, http=cached_http
            )
            cached_records = cached.collect(context(rare_variant()))

            unavailable_http = StubHttp(([],))
            unavailable_variant = rare_variant(
                change="A1V", spdi="NC_000007.14:10:A:G"
            )
            unavailable = AlphaMissenseEvidenceProvider(
                database=database, http=unavailable_http
            )
            self.assertEqual(unavailable.collect(context(unavailable_variant)), [])
            self.assertEqual(
                AlphaMissenseEvidenceProvider(
                    database=database, http=StubHttp()
                ).collect(context(unavailable_variant)),
                [],
            )

            with sqlite3.connect(database) as connection:
                rows = connection.execute(
                    "SELECT COUNT(*) FROM alphamissense_variant_cache"
                ).fetchone()[0]

        self.assertEqual(first_records, cached_records)
        self.assertEqual(len(first_http.posts), 1)
        self.assertEqual(cached_http.posts, [])
        self.assertEqual(len(unavailable_http.posts), 1)
        self.assertEqual(rows, 2)
        self.assertIn("1 variants from SQLite", cached.last_message or "")

    def test_protein_mismatch_is_not_attributed_and_is_cached(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "runs.sqlite3"
            http = StubHttp((vep_response(amino_acids="L/Q"),))
            provider = AlphaMissenseEvidenceProvider(
                database=database, http=http
            )
            self.assertEqual(provider.collect(context(rare_variant())), [])
            self.assertIn(
                "no exact protein match", provider.last_message or ""
            )
            second_http = StubHttp()
            second = AlphaMissenseEvidenceProvider(
                database=database, http=second_http
            )
            self.assertEqual(second.collect(context(rare_variant())), [])
        self.assertEqual(second_http.posts, [])

    def test_non_snv_is_skipped_without_web_request(self):
        http = StubHttp()
        provider = AlphaMissenseEvidenceProvider(http=http)
        deletion = rare_variant(
            change="E746_A750del",
            spdi="NC_000007.14:55190000:ACGT:A",
        )
        self.assertEqual(provider.collect(context(deletion)), [])
        self.assertEqual(http.posts, [])
        self.assertIn("not an unambiguous GRCh38 SNV", provider.last_message or "")


class AlphaMissenseExportTests(unittest.TestCase):
    def test_shared_and_evidence_csv_expose_prediction_columns(self):
        prediction = EvidenceRecord(
            source="alphamissense",
            source_record_id="variant|transcript",
            evidence_layer="computational_prediction",
            gene="EGFR",
            retrieved_at="2026-08-16T00:00:00+00:00",
            variant_key="GRCh38:7:55191822:T>G",
            transcript="ENST00000275493",
            protein_change="L858R",
            prediction_score=0.9968,
            prediction_class="likely_pathogenic",
        )
        shared = SharedMutation(
            gene="EGFR",
            protein_change="L858R",
            diseases=("Rare EGFR disorder",),
            orphanet_codes=("999",),
            clinvar_ids=("101",),
            clinvar_accessions=("VCV101",),
            germline_classifications=("Pathogenic",),
            canonical_spdis=("NC_000007.14:55191821:T:G",),
            tcga_projects=("TCGA-LUAD",),
            cancer_patient_count=1,
            cancer_event_count=1,
            genomic_match_count=1,
            genomic_confirmations=("GRCh38:7:55191822:T>G",),
            ssm_ids=("ssm-1",),
            match_level="confirmed",
        )
        result = AnalysisResult(
            gene="EGFR",
            generated_at="2026-08-16T00:00:00+00:00",
            gdc_program="TCGA",
            gdc_projects=(),
            min_cancer_patients=1,
            excluded_conditions=(),
            clinvar_record_limit=None,
            allowed_germline_classifications=("Pathogenic",),
            gdc_data_release=None,
            rare_variant_count=1,
            cancer_observation_count=1,
            shared_mutation_count=1,
            shared_mutations=(shared,),
            mutation_plot_data=(
                MutationPlotRecord(
                    gene="EGFR",
                    protein_change="L858R",
                    amino_acid_position=858,
                    in_rare_disease=True,
                    in_oncology=True,
                    diseases=("Rare EGFR disorder",),
                    orphanet_codes=("999",),
                    germline_classifications=("Pathogenic",),
                    cancer_patient_count=1,
                    cancer_event_count=1,
                    project_patient_counts=(),
                    match_level="confirmed",
                ),
            ),
            evidence_record_count=1,
            evidence_records=(prediction,),
        )
        with tempfile.TemporaryDirectory() as directory:
            shared_path = write_csv(result, Path(directory) / "shared.csv")
            evidence_path = write_evidence_csv(
                result, Path(directory) / "evidence.csv"
            )
            with shared_path.open(encoding="utf-8") as handle:
                shared_row = next(csv.DictReader(handle))
            with evidence_path.open(encoding="utf-8") as handle:
                evidence_row = next(csv.DictReader(handle))
        self.assertEqual(shared_row["alphamissense_scores"], "0.9968")
        self.assertEqual(
            shared_row["alphamissense_classes"], "likely_pathogenic"
        )
        self.assertEqual(evidence_row["prediction_score"], "0.9968")


if __name__ == "__main__":
    unittest.main()
