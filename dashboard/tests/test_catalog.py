from __future__ import annotations

import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

from rare_oncology_dashboard.catalog import (
    ANALYSIS_KIND,
    THERAPY_KIND,
    SQLiteCatalog,
    discover_databases,
)


SCHEMA = """
CREATE TABLE execution_cache (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    kind TEXT NOT NULL,
    gene TEXT NOT NULL,
    parameter_hash TEXT NOT NULL,
    parameters_json TEXT NOT NULL,
    result_json TEXT NOT NULL,
    created_at TEXT NOT NULL,
    last_accessed_at TEXT NOT NULL,
    cache_hit_count INTEGER NOT NULL,
    schema_version INTEGER NOT NULL
)
"""

STRING_SCHEMA = """
CREATE TABLE string_network_cache (
    gene TEXT NOT NULL,
    species INTEGER NOT NULL,
    string_version TEXT NOT NULL,
    required_score INTEGER NOT NULL,
    network_type TEXT NOT NULL,
    added_nodes INTEGER NOT NULL,
    result_json TEXT NOT NULL,
    created_at TEXT NOT NULL,
    last_accessed_at TEXT NOT NULL,
    cache_hit_count INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (
        gene, species, string_version, required_score, network_type, added_nodes
    )
)
"""


def make_database(path: Path, rows: list[tuple[str, str, str, dict]]) -> None:
    connection = sqlite3.connect(path)
    connection.execute(SCHEMA)
    for kind, gene, created_at, result in rows:
        connection.execute(
            """
            INSERT INTO execution_cache (
              kind, gene, parameter_hash, parameters_json, result_json,
              created_at, last_accessed_at, cache_hit_count, schema_version
            ) VALUES (?, ?, ?, ?, ?, ?, ?, 0, 1)
            """,
            (
                kind,
                gene,
                f"hash-{gene}-{created_at}",
                json.dumps({"request": gene}),
                json.dumps(result),
                created_at,
                created_at,
            ),
        )
    connection.commit()
    connection.close()


class SQLiteCatalogTests(unittest.TestCase):
    def test_discovers_and_federates_databases(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "a" / "rare_oncology.sqlite3"
            second = root / "b" / "rare_oncology.sqlite3"
            first.parent.mkdir()
            second.parent.mkdir()
            make_database(
                first,
                [(ANALYSIS_KIND, "EGFR", "2026-01-01T00:00:00+00:00", {"gene": "EGFR", "v": 1})],
            )
            make_database(
                second,
                [
                    (ANALYSIS_KIND, "EGFR", "2026-02-01T00:00:00+00:00", {"gene": "EGFR", "v": 2}),
                    (THERAPY_KIND, "EGFR", "2026-02-02T00:00:00+00:00", {"gene": "EGFR"}),
                ],
            )
            paths = discover_databases(root, root / "master.sqlite3")
            catalog = SQLiteCatalog(paths)

            self.assertEqual(paths, (first.resolve(), second.resolve()))
            self.assertEqual(len(catalog.genes()), 1)
            entry = catalog.genes()[0]
            self.assertEqual(entry.analysis_runs, 2)
            self.assertEqual(entry.therapy_runs, 1)
            latest = catalog.latest("egfr", ANALYSIS_KIND)
            self.assertIsNotNone(latest)
            assert latest is not None
            self.assertEqual(latest.result["v"], 2)
            self.assertEqual(latest.database_path, second.resolve())
            self.assertEqual(latest.schema_version, 1)

            all_analyses = catalog.all_snapshots(ANALYSIS_KIND)
            self.assertEqual([snapshot.result["v"] for snapshot in all_analyses], [2, 1])
            self.assertEqual(
                [snapshot.database_path for snapshot in all_analyses],
                [second.resolve(), first.resolve()],
            )

    def test_ignores_unrelated_sqlite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            unrelated = root / "other.sqlite3"
            connection = sqlite3.connect(unrelated)
            connection.execute("CREATE TABLE other (id INTEGER)")
            connection.close()
            catalog = SQLiteCatalog((unrelated,))
            self.assertEqual(catalog.available_databases(), ())
            self.assertIn("execution_cache", catalog.warnings()[0])

    def test_discovers_tracked_bootstrap_from_an_additional_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            results = root / "results"
            bootstrap = root / "bootstrap"
            results.mkdir()
            bootstrap.mkdir()
            seed = bootstrap / "rare_oncology_demo.sqlite3"
            make_database(
                seed,
                [(ANALYSIS_KIND, "PDGFRB", "2026-01-01T00:00:00+00:00", {"gene": "PDGFRB"})],
            )

            paths = discover_databases(
                results,
                results / "master.sqlite3",
                additional_roots=(bootstrap,),
            )

            self.assertEqual(paths, (seed.resolve(),))
            self.assertEqual(SQLiteCatalog(paths).genes()[0].gene, "PDGFRB")

    def test_reads_latest_string_network_even_from_string_only_database(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "string.sqlite3"
            connection = sqlite3.connect(database)
            connection.execute(STRING_SCHEMA)
            for score, created_at in (
                (400, "2026-01-01T00:00:00+00:00"),
                (700, "2026-02-01T00:00:00+00:00"),
            ):
                connection.execute(
                    """
                    INSERT INTO string_network_cache VALUES (
                      'EGFR', 9606, '12.0', ?, 'functional', 20, ?, ?, ?, 0
                    )
                    """,
                    (
                        score,
                        json.dumps({"gene": "EGFR", "required_score": score}),
                        created_at,
                        created_at,
                    ),
                )
            connection.commit()
            connection.close()

            catalog = SQLiteCatalog((database,))
            self.assertEqual(catalog.available_databases(), (database.resolve(),))
            self.assertEqual(catalog.genes(), ())
            snapshots = catalog.string_networks("egfr")
            self.assertEqual(len(snapshots), 2)
            self.assertEqual(snapshots[0].required_score, 700)
            self.assertEqual(catalog.latest_string_network("EGFR"), snapshots[0])


if __name__ == "__main__":
    unittest.main()
