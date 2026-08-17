from __future__ import annotations

import csv
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from rare_oncology.cli import build_parser, main
from rare_oncology.models import AnalysisResult
from rare_oncology.plots import (
    create_string_network_plot,
    generate_string_network_plots,
)
from rare_oncology.string_network import (
    STRING_SPECIES_HUMAN,
    STRING_VERSION,
    StringNetworkError,
    StringNetworkRun,
    get_string_network,
    write_string_edges_csv,
    write_string_network_json,
    write_string_nodes_csv,
)
from rare_oncology.string_cli import (
    build_parser as build_string_parser,
    main as string_main,
)


MAPPING = [
    {
        "queryItem": "EGFR",
        "queryIndex": 0,
        "stringId": "9606.ENSP00000275493",
        "ncbiTaxonId": 9606,
        "taxonName": "Homo sapiens",
        "preferredName": "EGFR",
        "annotation": "Receptor tyrosine-protein kinase erbB-1",
    }
]

NETWORK = [
    {
        "stringId_A": "9606.ENSP00000275493",
        "stringId_B": "9606.ENSP00000266970",
        "preferredName_A": "EGFR",
        "preferredName_B": "GRB2",
        "ncbiTaxonId": 9606,
        "score": 0.999,
        "nscore": 0.0,
        "fscore": 0.0,
        "pscore": 0.0,
        "ascore": 0.81,
        "escore": 0.94,
        "dscore": 0.9,
        "tscore": 0.88,
    },
    {
        "stringId_A": "9606.ENSP00000275493",
        "stringId_B": "9606.ENSP00000344456",
        "preferredName_A": "EGFR",
        "preferredName_B": "ERBB2",
        "score": 0.995,
        "nscore": 0.0,
        "fscore": 0.0,
        "pscore": 0.0,
        "ascore": 0.71,
        "escore": 0.91,
        "dscore": 0.9,
        "tscore": 0.78,
    },
    {
        "stringId_A": "9606.ENSP00000266970",
        "stringId_B": "9606.ENSP00000344456",
        "preferredName_A": "GRB2",
        "preferredName_B": "ERBB2",
        "score": 0.72,
        "nscore": 0.0,
        "fscore": 0.0,
        "pscore": 0.0,
        "ascore": 0.4,
        "escore": 0.2,
        "dscore": 0.5,
        "tscore": 0.6,
    },
]


class StubHttp:
    def __init__(self, responses=()) -> None:
        self.responses = iter(responses)
        self.gets: list[tuple[str, dict]] = []

    def get_json(self, url, params):
        self.gets.append((url, params))
        return next(self.responses)


def network_run(*, database_path: str | None = None) -> StringNetworkRun:
    return get_string_network(
        "EGFR",
        database=database_path,
        required_score=700,
        added_nodes=20,
        http=StubHttp((MAPPING, NETWORK)),
    )


def analysis_result() -> AnalysisResult:
    return AnalysisResult(
        gene="EGFR",
        generated_at="2026-08-16T00:00:00+00:00",
        gdc_program=None,
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
    )


class StringNetworkTests(unittest.TestCase):
    def test_focused_cli_exports_from_cache_path_without_analysis(self) -> None:
        parsed = build_string_parser().parse_args(
            ["TGFBR1", "--string-score", "800", "--string-partners", "30"]
        )
        self.assertEqual(parsed.gene, "TGFBR1")
        self.assertEqual(parsed.string_score, 800)
        self.assertEqual(parsed.string_partners, 30)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = network_run()
            with patch(
                "rare_oncology.string_cli.get_string_network", return_value=run
            ) as fetch:
                self.assertEqual(
                    string_main(
                        [
                            "EGFR",
                            "--database",
                            str(root / "master.sqlite3"),
                            "--output-dir",
                            str(root / "string"),
                            "--offline",
                        ]
                    ),
                    0,
                )
            self.assertTrue(fetch.call_args.kwargs["offline"])
            self.assertTrue((root / "string" / "egfr_string_network.json").is_file())
            self.assertTrue((root / "string" / "egfr_string_nodes.csv").is_file())
            self.assertTrue((root / "string" / "egfr_string_edges.csv").is_file())

    def test_maps_human_gene_then_fetches_versioned_network(self) -> None:
        http = StubHttp((MAPPING, NETWORK))
        run = get_string_network(
            "egfr",
            required_score=700,
            network_type="functional",
            added_nodes=20,
            http=http,
        )

        self.assertFalse(run.cache_hit)
        self.assertEqual(run.network.gene, "EGFR")
        self.assertEqual(run.network.species, STRING_SPECIES_HUMAN)
        self.assertEqual(run.network.string_version, STRING_VERSION)
        self.assertEqual(len(run.network.nodes), 3)
        self.assertEqual(len(run.network.edges), 3)
        self.assertTrue(run.network.nodes[0].is_query)
        self.assertIn("version-12-0.string-db.org", http.gets[0][0])
        self.assertEqual(http.gets[0][1]["species"], 9606)
        self.assertEqual(http.gets[1][1]["identifiers"], "9606.ENSP00000275493")
        self.assertEqual(http.gets[1][1]["required_score"], 700)
        self.assertEqual(http.gets[1][1]["add_nodes"], 20)

    def test_sqlite_is_consulted_before_web_and_tracks_hit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "runs.sqlite3"
            first_http = StubHttp((MAPPING, NETWORK))
            first = get_string_network("EGFR", database=database, http=first_http)
            cached_http = StubHttp()
            second = get_string_network(
                "EGFR", database=database, offline=True, http=cached_http
            )
            with sqlite3.connect(database) as connection:
                row = connection.execute(
                    "SELECT cache_hit_count, result_json FROM string_network_cache"
                ).fetchone()

        self.assertFalse(first.cache_hit)
        self.assertTrue(second.cache_hit)
        self.assertEqual(first.network, second.network)
        self.assertEqual(cached_http.gets, [])
        self.assertEqual(row[0], 1)
        self.assertEqual(json.loads(row[1])["gene"], "EGFR")

    def test_parameter_change_does_not_reuse_incompatible_network(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "runs.sqlite3"
            get_string_network(
                "EGFR",
                database=database,
                added_nodes=10,
                http=StubHttp((MAPPING, NETWORK)),
            )
            other = StubHttp((MAPPING, NETWORK))
            result = get_string_network(
                "EGFR", database=database, added_nodes=30, http=other
            )
        self.assertFalse(result.cache_hit)
        self.assertEqual(len(other.gets), 2)

    def test_offline_cache_miss_never_calls_web(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            http = StubHttp()
            with self.assertRaisesRegex(StringNetworkError, "SQLite"):
                get_string_network(
                    "EGFR",
                    database=Path(directory) / "runs.sqlite3",
                    offline=True,
                    http=http,
                )
        self.assertEqual(http.gets, [])

    def test_invalid_score_is_rejected_instead_of_serializing_nan(self) -> None:
        malformed = [dict(NETWORK[0], score=float("nan"))]
        with self.assertRaisesRegex(StringNetworkError, "out of range"):
            get_string_network("EGFR", http=StubHttp((MAPPING, malformed)))

    def test_exports_json_nodes_and_edges_with_safe_csv(self) -> None:
        mapping = [dict(MAPPING[0], annotation="=HYPERLINK(\"bad\")")]
        run = get_string_network("EGFR", http=StubHttp((mapping, NETWORK)))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            json_path = write_string_network_json(run.network, root / "network.json")
            nodes_path = write_string_nodes_csv(run.network, root / "nodes.csv")
            edges_path = write_string_edges_csv(run.network, root / "edges.csv")
            with nodes_path.open(encoding="utf-8") as handle:
                nodes = list(csv.DictReader(handle))

            self.assertEqual(json.loads(json_path.read_text())["gene"], "EGFR")
            self.assertEqual(len(nodes), 3)
            self.assertTrue(nodes[0]["annotation"].startswith("'="))
            with edges_path.open(encoding="utf-8") as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), 3)

    def test_plot_and_export_are_deterministic_and_nonempty(self) -> None:
        run = network_run()
        figure = create_string_network_plot(run.network)
        self.addCleanup(__import__("matplotlib.pyplot").pyplot.close, figure)
        self.assertIn("STRING v12.0", figure.axes[0].get_title())
        labels = {text.get_text() for text in figure.axes[0].texts}
        self.assertTrue({"EGFR", "GRB2", "ERBB2"} <= labels)
        with tempfile.TemporaryDirectory() as directory:
            paths = generate_string_network_plots(run.network, directory)
            self.assertEqual({path.suffix for path in paths}, {".png", ".svg"})
            self.assertTrue(all(path.stat().st_size > 100 for path in paths))

    def test_cli_exports_network_and_requests_plot_when_selected(self) -> None:
        parsed = build_parser().parse_args(
            [
                "EGFR",
                "--string-network",
                "--string-partners",
                "15",
                "--string-score",
                "800",
                "--string-network-type",
                "physical",
            ]
        )
        self.assertTrue(parsed.string_network)
        self.assertEqual(parsed.string_partners, 15)
        self.assertEqual(parsed.string_score, 800)
        self.assertEqual(parsed.string_network_type, "physical")

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "runs.sqlite3"
            output = root / "output"
            cached_analysis = SimpleNamespace(
                result=analysis_result(),
                cache_hit=True,
                execution_id=1,
                database_path=str(database),
            )
            run = network_run()
            with (
                patch(
                    "rare_oncology.cli.find_shared_mutations_cached",
                    return_value=cached_analysis,
                ),
                patch("rare_oncology.cli.get_string_network", return_value=run) as string_call,
                patch("rare_oncology.cli.generate_gene_plots", return_value={}),
                patch(
                    "rare_oncology.cli.generate_string_network_plots",
                    return_value=(output / "plots" / "egfr_string_network.png",),
                ) as plot_call,
            ):
                self.assertEqual(
                    main(
                        [
                            "EGFR",
                            "--database",
                            str(database),
                            "--output-dir",
                            str(output),
                            "--offline",
                            "--string-network",
                            "--plots",
                        ]
                    ),
                    0,
                )

            self.assertTrue((output / "egfr_string_network.json").is_file())
            self.assertTrue((output / "egfr_string_nodes.csv").is_file())
            self.assertTrue((output / "egfr_string_edges.csv").is_file())
            self.assertTrue(string_call.call_args.kwargs["offline"])
            plot_call.assert_called_once()


if __name__ == "__main__":
    unittest.main()
