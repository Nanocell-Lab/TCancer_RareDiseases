from __future__ import annotations

import json
import sqlite3
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

from rare_oncology.completion import (
    complete_missing_queries,
    discover_completion_databases,
)
from rare_oncology.completion_cli import main
from rare_oncology.evidence import ONCOLOGY_EVIDENCE
from rare_oncology.models import (
    AnalysisResult,
    EvidenceRecord,
    SourceStatus,
    TherapyAssociation,
)
from rare_oncology.protein_annotations import ProteinAnnotation
from rare_oncology.storage import ANALYSIS_KIND, ExecutionDatabase
from rare_oncology.therapy import VARIANT_THERAPY_EVIDENCE


def analysis_result() -> AnalysisResult:
    records = (
        EvidenceRecord(
            source="clinvar",
            source_record_id="VCV1|ORPHA:1",
            evidence_layer="variant_assertion",
            gene="EGFR",
            retrieved_at="2026-08-17T00:00:00+00:00",
            variant_key="GRCh38:7:55191822:T:G",
            assembly="GRCh38",
            chromosome="7",
            position=55191822,
            reference="T",
            alternate="G",
            protein_change="L858R",
            disease_id="ORPHA:1",
            disease_name="Rare EGFR disorder",
            classification="Pathogenic",
            details=(("clinvar_id", "1"),),
        ),
        EvidenceRecord(
            source="gdc",
            source_record_id="ssm-1|case-1",
            evidence_layer="case_observation",
            gene="EGFR",
            retrieved_at="2026-08-17T00:00:00+00:00",
            variant_key="GRCh38:7:55191822:T:G",
            assembly="GRCh38",
            chromosome="7",
            position=55191822,
            reference="T",
            alternate="G",
            protein_change="L858R",
            details=(
                ("case_id", "case-1"),
                ("patient_id", "patient-1"),
                ("project_id", "TCGA-LUAD"),
                ("ssm_id", "ssm-1"),
            ),
        ),
    )
    return AnalysisResult(
        gene="EGFR",
        generated_at="2026-08-17T00:00:00+00:00",
        gdc_program="TCGA",
        gdc_projects=("TCGA-LUAD",),
        min_cancer_patients=1,
        excluded_conditions=(),
        clinvar_record_limit=None,
        allowed_germline_classifications=("Pathogenic",),
        gdc_data_release="46",
        rare_variant_count=1,
        cancer_observation_count=1,
        shared_mutation_count=0,
        shared_mutations=(),
        mutation_plot_data=(),
        evidence_record_count=len(records),
        evidence_records=records,
        source_statuses=(
            SourceStatus("clinvar", "ok", "variant_assertion", 1),
            SourceStatus("gdc", "ok", "case_observation", 1),
        ),
    )


class StubEvidenceProvider:
    name = "mockevidence"
    evidence_layer = ONCOLOGY_EVIDENCE
    url = "https://example.test/evidence"
    source_version = "1"

    def __init__(self) -> None:
        self.calls = 0
        self.context_sizes: tuple[int, int] | None = None

    def cache_identity(self):
        return {"url": self.url, "version": self.source_version}

    def collect(self, context):
        self.calls += 1
        self.context_sizes = (len(context.rare_variants), len(context.cancer_observations))
        return [
            EvidenceRecord(
                source=self.name,
                source_record_id="mock-1",
                evidence_layer=self.evidence_layer,
                gene=context.gene,
                retrieved_at="2026-08-17T00:00:00+00:00",
                protein_change="L858R",
            )
        ]


class StubTherapyProvider:
    name = "mocktherapy"
    layer = VARIANT_THERAPY_EVIDENCE
    url = "https://example.test/therapy"
    source_version = "1"

    def __init__(self) -> None:
        self.calls = 0

    def cache_identity(self):
        return {"url": self.url, "version": self.source_version}

    def collect(self, context):
        self.calls += 1
        return [
            TherapyAssociation(
                layer=self.layer,
                source=self.name,
                source_record_id="therapy-1",
                gene=context.gene,
                association_scope="gene_target",
                therapy_name="Mock drug",
                therapy_components=(),
                therapy_status="experimental",
                retrieved_at="2026-08-17T00:00:00+00:00",
            )
        ]


class FailingEvidenceProvider(StubEvidenceProvider):
    name = "failingevidence"

    def collect(self, context):
        self.calls += 1
        raise ValueError("temporary source failure")


class StubStringHttp:
    def __init__(self) -> None:
        self.calls = 0

    def get_json(self, url, params):
        self.calls += 1
        if url.endswith("get_string_ids"):
            return [
                {
                    "stringId": "9606.ENSP00000275493",
                    "preferredName": "EGFR",
                    "annotation": "EGFR receptor",
                }
            ]
        return [
            {
                "stringId_A": "9606.ENSP00000275493",
                "stringId_B": "9606.ENSP00000354587",
                "preferredName_A": "EGFR",
                "preferredName_B": "ERBB2",
                "score": 0.9,
                "nscore": 0.0,
                "fscore": 0.0,
                "pscore": 0.0,
                "ascore": 0.1,
                "escore": 0.8,
                "dscore": 0.8,
                "tscore": 0.4,
            }
        ]


class StubProteinClient:
    def __init__(self) -> None:
        self.calls = 0

    def fetch(self, gene):
        self.calls += 1
        return ProteinAnnotation(
            gene=gene,
            accession="P00533",
            length=4,
            domains=(),
            secondary_structure=(),
            source_url="https://rest.uniprot.org/uniprotkb/P00533.json",
            sequence="ACDE",
        )


class CompletionTests(unittest.TestCase):
    def test_backfill_queries_each_missing_layer_once_then_reuses_sqlite(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "runs.sqlite3"
            with ExecutionDatabase(database) as db:
                db.save(ANALYSIS_KIND, "EGFR", {}, analysis_result())
            evidence = StubEvidenceProvider()
            therapy = StubTherapyProvider()
            string_http = StubStringHttp()
            protein = StubProteinClient()
            kwargs = {
                "sources": ("mockevidence", "mocktherapy", "string", "uniprot"),
                "evidence_providers": (evidence,),
                "therapy_providers": (therapy,),
                "string_http": string_http,
                "protein_client": protein,
            }
            first = complete_missing_queries((database,), **kwargs)
            second = complete_missing_queries((database,), **kwargs)
            with ExecutionDatabase(database) as db:
                analyses = db.list_executions(gene="EGFR", kind=ANALYSIS_KIND)
                therapies = db.list_executions(gene="EGFR", kind="therapy_analysis")

        self.assertEqual(first.query_count, 4)
        self.assertEqual(first.error_count, 0)
        self.assertEqual(second.query_count, 0)
        self.assertEqual(second.cache_hit_count, 4)
        self.assertEqual(evidence.calls, 1)
        self.assertEqual(evidence.context_sizes, (1, 1))
        self.assertEqual(therapy.calls, 1)
        self.assertEqual(string_http.calls, 2)
        self.assertEqual(protein.calls, 1)
        self.assertEqual(len(analyses), 2)
        self.assertEqual(analyses[0].parameter_hash, analyses[1].parameter_hash)
        self.assertEqual(len(therapies), 1)

    def test_offline_and_dry_run_never_query_missing_source(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "runs.sqlite3"
            with ExecutionDatabase(database) as db:
                db.save(ANALYSIS_KIND, "EGFR", {}, analysis_result())
            http = StubStringHttp()
            offline = complete_missing_queries(
                (database,), sources=("string",), offline=True, string_http=http
            )
            dry = complete_missing_queries(
                (database,), sources=("string", "alphamissense"), dry_run=True,
                string_http=http,
            )
            with sqlite3.connect(database) as connection:
                completion_table = connection.execute(
                    "SELECT 1 FROM sqlite_master WHERE name='completion_query_status'"
                ).fetchone()
                alphamissense_table = connection.execute(
                    "SELECT 1 FROM sqlite_master WHERE name='alphamissense_variant_cache'"
                ).fetchone()
        self.assertEqual(http.calls, 0)
        self.assertEqual(offline.missing_count, 1)
        self.assertEqual(dry.missing_count, 2)
        self.assertEqual(dry.items[0].action, "would_query")
        self.assertIsNone(completion_table)
        self.assertIsNone(alphamissense_table)

    def test_errors_are_audited_and_only_retried_explicitly(self):
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "runs.sqlite3"
            with ExecutionDatabase(database) as db:
                db.save(ANALYSIS_KIND, "EGFR", {}, analysis_result())
            provider = FailingEvidenceProvider()
            kwargs = {
                "sources": (provider.name,),
                "evidence_providers": (provider,),
            }
            first = complete_missing_queries((database,), **kwargs)
            second = complete_missing_queries((database,), **kwargs)
            third = complete_missing_queries(
                (database,), retry_errors=True, **kwargs
            )
        self.assertEqual(first.error_count, 1)
        self.assertEqual(second.items[0].action, "cache_hit")
        self.assertEqual(third.items[0].action, "queried")
        self.assertEqual(provider.calls, 2)

    def test_discovery_and_cli_dry_run_write_machine_readable_report(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "nested" / "runs.sqlite3"
            database.parent.mkdir()
            with ExecutionDatabase(database) as db:
                db.save(ANALYSIS_KIND, "EGFR", {}, analysis_result())
            self.assertEqual(discover_completion_databases(root), (database.resolve(),))
            output = root / "completion"
            with redirect_stdout(StringIO()):
                code = main(
                    [
                        "--results-root", str(root), "--source", "string",
                        "--dry-run", "--output-dir", str(output),
                    ]
                )
            payload = json.loads(
                (output / "completion_report.json").read_text(encoding="utf-8")
            )
        self.assertEqual(code, 0)
        self.assertEqual(payload["gene_count"], 1)
        self.assertEqual(payload["missing_count"], 1)


if __name__ == "__main__":
    unittest.main()
