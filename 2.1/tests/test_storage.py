from __future__ import annotations

import tempfile
import unittest
import sqlite3
import os
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime, timedelta, timezone
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from rare_oncology.models import (
    AnalysisResult,
    EvidenceRecord,
    MutationPlotRecord,
    SharedMutation,
    SourceStatus,
    TherapyAssociation,
    TherapyResult,
)
from rare_oncology.cli import build_parser, main
from rare_oncology.storage import (
    ANALYSIS_KIND,
    CacheMissError,
    ExecutionDatabase,
    analysis_result_from_dict,
    find_shared_mutations_cached,
    find_therapy_associations_cached,
)
from rare_oncology.therapy import (
    OncoKBTherapyProvider,
    TherapyContext,
    VARIANT_THERAPY_EVIDENCE,
)


def analysis_result(change: str = "L858R") -> AnalysisResult:
    shared = SharedMutation(
        gene="EGFR",
        protein_change=change,
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
        other_inherited_diseases=("Rare EGFR disorder",),
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
        gdc_data_release="46",
        rare_variant_count=1,
        cancer_observation_count=1,
        shared_mutation_count=1,
        shared_mutations=(shared,),
        mutation_plot_data=(
            MutationPlotRecord(
                gene="EGFR",
                protein_change=change,
                amino_acid_position=858,
                in_rare_disease=True,
                in_oncology=True,
                diseases=("Rare EGFR disorder",),
                orphanet_codes=("999",),
                germline_classifications=("Pathogenic",),
                cancer_patient_count=1,
                cancer_event_count=1,
                project_patient_counts=(("TCGA-LUAD", 1),),
                match_level="confirmed",
                disease_variant_counts=(("Rare EGFR disorder", 1),),
                other_inherited_diseases=("Rare EGFR disorder",),
            ),
        ),
        evidence_record_count=1,
        evidence_records=(
            EvidenceRecord(
                source="fixture",
                source_record_id="1",
                evidence_layer="variant_assertion",
                gene="EGFR",
                retrieved_at="2026-08-11T00:00:00+00:00",
                details=(("key", "value"),),
            ),
        ),
        source_statuses=(
            SourceStatus(
                source="fixture",
                status="ok",
                evidence_layer="variant_assertion",
                record_count=1,
            ),
        ),
    )


class MutableClock:
    def __init__(self) -> None:
        self.value = datetime(2026, 8, 11, tzinfo=timezone.utc)

    def __call__(self) -> datetime:
        return self.value


class DatabaseTests(unittest.TestCase):
    def test_legacy_snapshot_rehydrates_with_separate_disease_categories(self):
        payload = analysis_result().to_dict()
        shared = payload["shared_mutations"][0]
        plotted = payload["mutation_plot_data"][0]
        shared.pop("hereditary_cancer_diseases")
        shared.pop("other_inherited_diseases")
        plotted.pop("hereditary_cancer_diseases")
        plotted.pop("other_inherited_diseases")

        restored = analysis_result_from_dict(payload)

        self.assertEqual(
            restored.shared_mutations[0].other_inherited_diseases,
            ("Rare EGFR disorder",),
        )
        self.assertEqual(
            restored.mutation_plot_data[0].other_inherited_diseases,
            ("Rare EGFR disorder",),
        )

    def test_round_trip_preserves_nested_dataclasses_and_tracks_hits(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            with ExecutionDatabase(path) as database:
                execution_id = database.save(
                    ANALYSIS_KIND, "EGFR", {"program": "TCGA"}, analysis_result()
                )
                stored = database.load_latest(
                    ANALYSIS_KIND, "egfr", {"program": "TCGA"}
                )
                summaries = database.list_executions(gene="EGFR")
        self.assertIsNotNone(stored)
        assert stored is not None
        self.assertEqual(stored.execution_id, execution_id)
        self.assertIsInstance(stored.result, AnalysisResult)
        if not isinstance(stored.result, AnalysisResult):
            self.fail("Se esperaba AnalysisResult")
        analysis = stored.result
        self.assertEqual(analysis.shared_mutations[0].protein_change, "L858R")
        self.assertIsInstance(analysis.shared_mutations, tuple)
        self.assertEqual(
            analysis.mutation_plot_data[0].project_patient_counts,
            (("TCGA-LUAD", 1),),
        )
        self.assertEqual(
            analysis.mutation_plot_data[0].disease_variant_counts,
            (("Rare EGFR disorder", 1),),
        )
        self.assertEqual(
            analysis.shared_mutations[0].other_inherited_diseases,
            ("Rare EGFR disorder",),
        )
        self.assertEqual(
            analysis.mutation_plot_data[0].other_inherited_diseases,
            ("Rare EGFR disorder",),
        )
        self.assertEqual(analysis.evidence_records[0].details, (("key", "value"),))
        self.assertEqual(summaries[0].cache_hit_count, 1)

    def test_parameter_change_is_a_cache_miss(self):
        with ExecutionDatabase(":memory:") as database:
            database.save(ANALYSIS_KIND, "EGFR", {"program": "TCGA"}, analysis_result())
            stored = database.load_latest(
                ANALYSIS_KIND, "EGFR", {"program": "ALL"}
            )
        self.assertIsNone(stored)

    def test_expired_snapshot_is_not_returned(self):
        clock = MutableClock()
        with ExecutionDatabase(":memory:", clock=clock) as database:
            database.save(ANALYSIS_KIND, "EGFR", {}, analysis_result())
            clock.value += timedelta(hours=25)
            stored = database.load_latest(
                ANALYSIS_KIND, "EGFR", {}, max_age_hours=24
            )
        self.assertIsNone(stored)

    def test_nonfinite_ttl_is_rejected_by_api_and_cli(self):
        with ExecutionDatabase(":memory:") as database:
            for value in (float("nan"), float("inf")):
                with self.subTest(value=value), self.assertRaises(ValueError):
                    database.load_latest(
                        ANALYSIS_KIND, "EGFR", {}, max_age_hours=value
                    )
        for value in ("nan", "inf"):
            with self.subTest(cli=value), redirect_stderr(StringIO()), self.assertRaises(
                SystemExit
            ):
                build_parser().parse_args(
                    ["EGFR", "--cache-max-age-hours", value]
                )

    def test_sensitive_cache_parameters_are_not_persisted(self):
        class Provider:
            url = "https://example.org/provider?token=provider-secret"

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            with ExecutionDatabase(path) as database:
                database.save(
                    ANALYSIS_KIND,
                    "EGFR",
                    {
                        "api_key": "secret-key",
                        "email": "person@example.org",
                        "url": (
                            "https://alice:BasicSecret@example.org/data?"
                            "client_secret=QuerySecret&X-Amz-Signature=SigSecret"
                        ),
                        "provider": Provider(),
                    },
                    analysis_result(),
                )
            with sqlite3.connect(path) as connection:
                parameters = connection.execute(
                    "SELECT parameters_json FROM execution_cache"
                ).fetchone()[0]
        self.assertNotIn("secret-key", parameters)
        self.assertNotIn("person@example.org", parameters)
        self.assertNotIn("BasicSecret", parameters)
        self.assertNotIn("QuerySecret", parameters)
        self.assertNotIn("SigSecret", parameters)
        self.assertNotIn("provider-secret", parameters)

    @unittest.skipUnless(os.name == "posix", "permisos POSIX")
    def test_database_and_existing_sidecars_are_owner_only(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            with ExecutionDatabase(path) as database:
                database.save(ANALYSIS_KIND, "EGFR", {}, analysis_result())
                candidates = (path, Path(f"{path}-wal"), Path(f"{path}-shm"))
                existing = [candidate for candidate in candidates if candidate.exists()]
                self.assertIn(path, existing)
                for candidate in existing:
                    self.assertEqual(candidate.stat().st_mode & 0o777, 0o600)

    def test_file_content_change_invalidates_even_with_same_size_and_mtime(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "evidence.csv"
            path.write_text("AAAA", encoding="utf-8")
            original = path.stat()
            with ExecutionDatabase(":memory:") as database:
                database.save(
                    ANALYSIS_KIND,
                    "EGFR",
                    {"evidence_file": path},
                    analysis_result(),
                )
                path.write_text("BBBB", encoding="utf-8")
                path.touch()
                os.utime(path, ns=(original.st_atime_ns, original.st_mtime_ns))
                stored = database.load_latest(
                    ANALYSIS_KIND, "EGFR", {"evidence_file": path}
                )
        self.assertIsNone(stored)

    def test_save_rejects_kind_and_gene_mismatches(self):
        therapy = TherapyResult(
            gene="EGFR",
            generated_at="2026-08-11T00:00:00+00:00",
            protein_changes=(),
            variant_therapy_evidence=(),
            regulatory_therapy=(),
            trial_and_target_discovery=(),
            source_statuses=(),
        )
        with ExecutionDatabase(":memory:") as database:
            with self.assertRaises(TypeError):
                database.save(ANALYSIS_KIND, "EGFR", {}, therapy)
            with self.assertRaises(ValueError):
                database.save(ANALYSIS_KIND, "ALK", {}, analysis_result())
            self.assertEqual(database.list_executions(), ())


class CachedFunctionTests(unittest.TestCase):
    def test_analysis_queries_web_once_then_uses_database(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            with patch(
                "rare_oncology.storage.find_shared_mutations",
                return_value=analysis_result(),
            ) as web:
                first = find_shared_mutations_cached(
                    "EGFR",
                    database=path,
                    cache_parameters={"program": "TCGA"},
                )
                second = find_shared_mutations_cached(
                    "EGFR",
                    database=path,
                    cache_parameters={"program": "TCGA"},
                )
        self.assertFalse(first.cache_hit)
        self.assertTrue(second.cache_hit)
        self.assertEqual(first.execution_id, second.execution_id)
        self.assertEqual(web.call_count, 1)

    def test_refresh_executes_again_and_keeps_history(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            with patch(
                "rare_oncology.storage.find_shared_mutations",
                return_value=analysis_result(),
            ) as web:
                first = find_shared_mutations_cached(
                    "EGFR", database=path, cache_parameters={}
                )
                second = find_shared_mutations_cached(
                    "EGFR", database=path, cache_parameters={}, refresh=True
                )
            with ExecutionDatabase(path) as database:
                history = database.list_executions(kind=ANALYSIS_KIND)
        self.assertEqual(web.call_count, 2)
        self.assertNotEqual(first.execution_id, second.execution_id)
        self.assertEqual(len(history), 2)

    def test_explicit_cache_metadata_does_not_hide_argument_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            with patch(
                "rare_oncology.storage.find_shared_mutations",
                return_value=analysis_result(),
            ) as web:
                first = find_shared_mutations_cached(
                    "EGFR",
                    database=path,
                    cache_parameters={"label": "same"},
                    program="TCGA",
                )
                second = find_shared_mutations_cached(
                    "EGFR",
                    database=path,
                    cache_parameters={"label": "same"},
                    program=None,
                )
        self.assertFalse(first.cache_hit)
        self.assertFalse(second.cache_hit)
        self.assertEqual(web.call_count, 2)

    def test_biomedical_variant_key_is_not_redacted_or_collapsed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            with patch(
                "rare_oncology.storage.find_shared_mutations",
                return_value=analysis_result(),
            ) as web:
                first = find_shared_mutations_cached(
                    "EGFR",
                    database=path,
                    cache_parameters={"variant_key": "L858R"},
                )
                second = find_shared_mutations_cached(
                    "EGFR",
                    database=path,
                    cache_parameters={"variant_key": "T790M"},
                )
        self.assertFalse(first.cache_hit)
        self.assertFalse(second.cache_hit)
        self.assertEqual(web.call_count, 2)

    def test_one_shot_analysis_iterables_are_materialized_for_key_and_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            with patch(
                "rare_oncology.storage.find_shared_mutations",
                return_value=analysis_result(),
            ) as web:
                first = find_shared_mutations_cached(
                    "EGFR",
                    database=path,
                    projects=(item for item in ("TCGA-LUAD",)),
                )
                second = find_shared_mutations_cached(
                    "EGFR",
                    database=path,
                    projects=(item for item in ("TARGET-AML",)),
                )
        self.assertFalse(first.cache_hit)
        self.assertFalse(second.cache_hit)
        self.assertEqual(web.call_count, 2)
        self.assertEqual(web.call_args_list[0].kwargs["projects"], ("TCGA-LUAD",))
        self.assertEqual(web.call_args_list[1].kwargs["projects"], ("TARGET-AML",))

    def test_repeatable_nonsequence_iterables_have_distinct_cache_keys(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            with patch(
                "rare_oncology.storage.find_shared_mutations",
                return_value=analysis_result(),
            ) as web:
                first = find_shared_mutations_cached(
                    "EGFR",
                    database=path,
                    projects={"TCGA-LUAD": 1}.keys(),
                )
                second = find_shared_mutations_cached(
                    "EGFR",
                    database=path,
                    projects={"TARGET-AML": 1}.keys(),
                )
        self.assertFalse(first.cache_hit)
        self.assertFalse(second.cache_hit)
        self.assertEqual(web.call_count, 2)

    def test_evidence_provider_identity_distinguishes_custom_datasets(self):
        class Provider:
            name = "external"
            evidence_layer = "variant_assertion"
            url = "https://example.org/evidence"
            source_version = None

            def __init__(self, dataset: str) -> None:
                self.dataset = dataset

            def cache_identity(self):
                return {"url": self.url, "dataset": self.dataset}

            def collect(self, context):
                return []

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            with patch(
                "rare_oncology.storage.find_shared_mutations",
                side_effect=(analysis_result("L858R"), analysis_result("T790M")),
            ) as web:
                first = find_shared_mutations_cached(
                    "EGFR",
                    database=path,
                    evidence_providers=[Provider("A")],
                )
                second = find_shared_mutations_cached(
                    "EGFR",
                    database=path,
                    evidence_providers=[Provider("B")],
                )
        self.assertFalse(first.cache_hit)
        self.assertFalse(second.cache_hit)
        self.assertEqual(second.result.shared_mutations[0].protein_change, "T790M")
        self.assertEqual(web.call_count, 2)

    def test_cached_analysis_rejects_provider_without_explicit_identity(self):
        class Provider:
            name = "external"
            evidence_layer = "variant_assertion"
            source_version = None
            url = "https://example.org/evidence"

            def collect(self, context):
                return []

        with patch("rare_oncology.storage.find_shared_mutations") as web:
            with self.assertRaisesRegex(TypeError, "cache_identity"):
                find_shared_mutations_cached(
                    "EGFR",
                    database=":memory:",
                    evidence_providers=[Provider()],
                )
        web.assert_not_called()

    def test_offline_cache_miss_never_calls_web(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            with patch("rare_oncology.storage.find_shared_mutations") as web:
                with self.assertRaises(CacheMissError):
                    find_shared_mutations_cached(
                        "EGFR",
                        database=path,
                        cache_parameters={},
                        offline=True,
                    )
        web.assert_not_called()

    def test_therapy_provider_runs_once_then_uses_database(self):
        class Provider:
            name = "fixture"
            layer = VARIANT_THERAPY_EVIDENCE
            source_version = "1"
            url = None

            def __init__(self) -> None:
                self.calls = 0

            def cache_identity(self):
                return {"url": self.url, "source_version": self.source_version}

            def collect(self, context):
                self.calls += 1
                return [
                    TherapyAssociation(
                        layer=self.layer,
                        source=self.name,
                        source_record_id="E1",
                        gene=context.gene,
                        association_scope="candidate_protein_variant",
                        therapy_name="Cachetinib",
                        therapy_components=("Cachetinib",),
                        therapy_status="curated_predictive_evidence",
                        retrieved_at="2026-08-11T00:00:00+00:00",
                        therapy_type="small_molecule",
                        chembl_ids=("CHEMBL42",),
                        canonical_smiles="CCN",
                    )
                ]

        provider = Provider()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            first = find_therapy_associations_cached(
                analysis_result(),
                database=path,
                providers=[provider],
                cache_parameters={"providers": ["fixture"]},
            )
            second = find_therapy_associations_cached(
                analysis_result(),
                database=path,
                providers=[provider],
                cache_parameters={"providers": ["fixture"]},
                offline=True,
            )
        self.assertFalse(first.cache_hit)
        self.assertTrue(second.cache_hit)
        self.assertEqual(provider.calls, 1)
        self.assertEqual(
            second.result.variant_therapy_evidence[0].therapy_name,
            "Cachetinib",
        )
        cached = second.result.variant_therapy_evidence[0]
        self.assertEqual(cached.therapy_type, "small_molecule")
        self.assertEqual(cached.chembl_ids, ("CHEMBL42",))
        self.assertEqual(cached.canonical_smiles, "CCN")

    def test_oncokb_token_is_not_persisted_in_sqlite(self):
        class Http:
            def get_json(self, url, params, *, headers):
                return {"dataVersion": "v7.4", "treatments": []}

        token = "DIRECT-ONCOKB-SECRET"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "runs.sqlite3"
            find_therapy_associations_cached(
                analysis_result(),
                database=path,
                providers=[OncoKBTherapyProvider(token, http=Http())],
            )
            with sqlite3.connect(path) as connection:
                persisted = "\n".join(
                    str(value)
                    for row in connection.execute(
                        "SELECT parameters_json, result_json FROM execution_cache"
                    )
                    for value in row
                )
        self.assertNotIn(token, persisted)

    def test_therapy_cache_key_includes_molecular_signature(self):
        class Provider:
            name = "empty"
            layer = VARIANT_THERAPY_EVIDENCE
            source_version = None
            url = None

            def __init__(self) -> None:
                self.calls = 0

            def cache_identity(self):
                return {"name": self.name}

            def collect(self, context: TherapyContext) -> list[TherapyAssociation]:
                self.calls += 1
                return []

        provider = Provider()
        with ExecutionDatabase(":memory:") as database:
            find_therapy_associations_cached(
                analysis_result("L858R"), database=database, providers=[provider]
            )
            second = find_therapy_associations_cached(
                analysis_result("T790M"), database=database, providers=[provider]
            )
        self.assertFalse(second.cache_hit)
        self.assertEqual(provider.calls, 2)

    def test_therapy_metadata_does_not_hide_provider_configuration(self):
        class Provider:
            name = "configured"
            layer = VARIANT_THERAPY_EVIDENCE
            source_version = "1"
            url = None

            def __init__(self, max_studies: int) -> None:
                self.max_studies = max_studies
                self.calls = 0

            def cache_identity(self):
                return {"max_studies": self.max_studies}

            def collect(self, context: TherapyContext) -> list[TherapyAssociation]:
                self.calls += 1
                return []

        first_provider = Provider(10)
        second_provider = Provider(20)
        with ExecutionDatabase(":memory:") as database:
            first = find_therapy_associations_cached(
                analysis_result(),
                database=database,
                providers=[first_provider],
                cache_parameters={"label": "same"},
            )
            second = find_therapy_associations_cached(
                analysis_result(),
                database=database,
                providers=[second_provider],
                cache_parameters={"label": "same"},
            )
        self.assertFalse(first.cache_hit)
        self.assertFalse(second.cache_hit)
        self.assertEqual(first_provider.calls, 1)
        self.assertEqual(second_provider.calls, 1)

    def test_mutating_provider_configuration_causes_cache_miss(self):
        class Provider:
            name = "configured"
            layer = VARIANT_THERAPY_EVIDENCE
            source_version = None
            url = "https://example.org/api"

            def __init__(self) -> None:
                self.max_studies = 10
                self.calls = 0

            def cache_identity(self):
                return {"url": self.url, "max_studies": self.max_studies}

            def collect(self, context: TherapyContext) -> list[TherapyAssociation]:
                self.calls += 1
                return []

        provider = Provider()
        with ExecutionDatabase(":memory:") as database:
            first = find_therapy_associations_cached(
                analysis_result(), database=database, providers=[provider]
            )
            provider.max_studies = 20
            second = find_therapy_associations_cached(
                analysis_result(), database=database, providers=[provider]
            )
        self.assertFalse(first.cache_hit)
        self.assertFalse(second.cache_hit)
        self.assertEqual(provider.calls, 2)

    def test_mutable_provider_response_metadata_does_not_change_cache_identity(self):
        class Provider:
            name = "mutable"
            layer = VARIANT_THERAPY_EVIDENCE
            url = "https://example.org/current"

            def __init__(self) -> None:
                self.source_version = None
                self.calls = 0
                self._configured_url = self.url

            def cache_identity(self):
                return {"url": self._configured_url}

            def collect(self, context: TherapyContext) -> list[TherapyAssociation]:
                self.calls += 1
                self.source_version = "release-2026-08"
                self.url = "https://example.org/releases/2026-08"
                return []

        provider = Provider()
        with ExecutionDatabase(":memory:") as database:
            first = find_therapy_associations_cached(
                analysis_result(), database=database, providers=[provider]
            )
            second = find_therapy_associations_cached(
                analysis_result(),
                database=database,
                providers=[provider],
                offline=True,
            )
        self.assertFalse(first.cache_hit)
        self.assertTrue(second.cache_hit)
        self.assertEqual(provider.calls, 1)

    def test_distinct_initial_provider_endpoints_do_not_share_cache(self):
        class Provider:
            name = "custom"
            layer = VARIANT_THERAPY_EVIDENCE
            source_version = None

            def __init__(self, url: str, therapy_name: str) -> None:
                self.url = url
                self.therapy_name = therapy_name
                self.calls = 0

            def cache_identity(self):
                return {"url": self.url}

            def collect(self, context: TherapyContext) -> list[TherapyAssociation]:
                self.calls += 1
                return [
                    TherapyAssociation(
                        layer=self.layer,
                        source=self.name,
                        source_record_id=self.therapy_name,
                        gene=context.gene,
                        association_scope="candidate_protein_variant",
                        therapy_name=self.therapy_name,
                        therapy_components=(self.therapy_name,),
                        therapy_status="curated_predictive_evidence",
                        retrieved_at="2026-08-11T00:00:00+00:00",
                    )
                ]

        first_provider = Provider("https://a.example/api", "Therapy A")
        second_provider = Provider("https://b.example/api", "Therapy B")
        with ExecutionDatabase(":memory:") as database:
            first = find_therapy_associations_cached(
                analysis_result(), database=database, providers=[first_provider]
            )
            second = find_therapy_associations_cached(
                analysis_result(), database=database, providers=[second_provider]
            )
        self.assertFalse(first.cache_hit)
        self.assertFalse(second.cache_hit)
        self.assertEqual(first_provider.calls, 1)
        self.assertEqual(second_provider.calls, 1)
        self.assertEqual(
            second.result.variant_therapy_evidence[0].therapy_name, "Therapy B"
        )

    def test_cli_reuses_sqlite_before_web_and_supports_offline(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "results"
            database = Path(directory) / "runs.sqlite3"
            argv = [
                "EGFR",
                "--database",
                str(database),
                "--output-dir",
                str(output),
            ]
            with patch(
                "rare_oncology.storage.find_shared_mutations",
                return_value=analysis_result(),
            ) as web, patch(
                "rare_oncology.cli.generate_gene_plots", return_value={}
            ) as plots, redirect_stdout(StringIO()) as stdout:
                self.assertEqual(main(argv), 0)
                self.assertEqual(main(argv), 0)
                self.assertEqual(main([*argv, "--offline", "--plots"]), 0)
        self.assertEqual(web.call_count, 1)
        self.assertIsNotNone(plots.call_args.kwargs["protein_annotation"])
        self.assertIn("SQLite molecular: HIT", stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
