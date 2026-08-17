"""Read-only federation over rare-oncology 2.1 execution databases."""

from __future__ import annotations

import json
import sqlite3
from contextlib import closing
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping


ANALYSIS_KIND = "mutation_analysis"
THERAPY_KIND = "therapy_analysis"
KNOWN_KINDS = (ANALYSIS_KIND, THERAPY_KIND)


class CatalogError(RuntimeError):
    """Raised when a requested snapshot cannot be read safely."""


@dataclass(frozen=True)
class ExecutionSnapshot:
    database_path: Path
    execution_id: int
    kind: str
    gene: str
    parameter_hash: str
    parameters: Mapping[str, Any]
    result: Mapping[str, Any]
    created_at: str
    last_accessed_at: str
    cache_hit_count: int

    @property
    def identity(self) -> str:
        """Role: Return the stable identity of this stored record."""
        return f"{self.database_path} · #{self.execution_id} · {self.created_at}"


@dataclass(frozen=True)
class StringNetworkSnapshot:
    database_path: Path
    gene: str
    species: int
    string_version: str
    required_score: int
    network_type: str
    added_nodes: int
    result: Mapping[str, Any]
    created_at: str
    last_accessed_at: str
    cache_hit_count: int

    @property
    def identity(self) -> str:
        """Role: Return the stable identity of this stored record."""
        return (
            f"{self.database_path} · STRING v{self.string_version} · "
            f"score {self.required_score} · {self.added_nodes} partners"
        )


@dataclass(frozen=True)
class GeneCatalogEntry:
    gene: str
    analysis_runs: int
    therapy_runs: int
    latest_analysis_at: str | None
    latest_therapy_at: str | None
    databases: tuple[Path, ...]

    @property
    def has_analysis(self) -> bool:
        """Role: Report whether the catalog entry contains analysis snapshots."""
        return self.analysis_runs > 0

    @property
    def has_therapy(self) -> bool:
        """Role: Report whether the catalog entry contains therapy snapshots."""
        return self.therapy_runs > 0


def discover_databases(results_root: str | Path, master: str | Path) -> tuple[Path, ...]:
    """Return execution SQLite files under 2.1/results plus the master if present."""

    root = Path(results_root).expanduser().resolve()
    master_path = Path(master).expanduser().resolve()
    paths = {
        path.resolve()
        for path in root.rglob("*.sqlite3")
        if path.is_file() and not path.name.endswith(("-wal", "-shm"))
    }
    if master_path.is_file():
        paths.add(master_path)
    return tuple(sorted(paths, key=lambda path: str(path)))


class SQLiteCatalog:
    """Federate immutable execution snapshots without copying their payloads."""

    def __init__(self, database_paths: tuple[str | Path, ...]) -> None:
        """Role: Compute init for the catalog workflow."""
        self.database_paths = tuple(
            dict.fromkeys(Path(path).expanduser().resolve() for path in database_paths)
        )

    def available_databases(self) -> tuple[Path, ...]:
        """Role: Return the databases that contribute to this catalog entry."""
        return tuple(path for path in self.database_paths if _is_supported_database(path))

    def warnings(self) -> tuple[str, ...]:
        """Role: Return non-fatal catalog warnings for display."""
        warnings: list[str] = []
        for path in self.database_paths:
            if not path.exists():
                warnings.append(f"Database not found: {path}")
                continue
            try:
                if not _is_supported_database(path):
                    warnings.append(
                        "Ignored database without execution_cache or "
                        f"string_network_cache: {path}"
                    )
            except sqlite3.Error as exc:
                warnings.append(f"Could not inspect {path}: {exc}")
        return tuple(warnings)

    def genes(self) -> tuple[GeneCatalogEntry, ...]:
        """Role: Build the current gene catalog across all configured databases."""
        aggregate: dict[str, dict[str, Any]] = {}
        for path in self.available_databases():
            if not _has_table(path, "execution_cache"):
                continue
            with closing(_connect_read_only(path)) as connection:
                rows = connection.execute(
                    """
                    SELECT gene, kind, COUNT(*) AS run_count, MAX(created_at) AS latest
                    FROM execution_cache
                    WHERE kind IN (?, ?)
                    GROUP BY gene, kind
                    """,
                    KNOWN_KINDS,
                ).fetchall()
            for row in rows:
                gene = str(row["gene"]).strip().upper()
                values = aggregate.setdefault(
                    gene,
                    {
                        "analysis_runs": 0,
                        "therapy_runs": 0,
                        "latest_analysis_at": None,
                        "latest_therapy_at": None,
                        "databases": set(),
                    },
                )
                run_count = int(row["run_count"])
                latest = str(row["latest"])
                if str(row["kind"]) == ANALYSIS_KIND:
                    values["analysis_runs"] += run_count
                    values["latest_analysis_at"] = _newest_iso(
                        values["latest_analysis_at"], latest
                    )
                else:
                    values["therapy_runs"] += run_count
                    values["latest_therapy_at"] = _newest_iso(
                        values["latest_therapy_at"], latest
                    )
                values["databases"].add(path)
        return tuple(
            GeneCatalogEntry(
                gene=gene,
                analysis_runs=int(values["analysis_runs"]),
                therapy_runs=int(values["therapy_runs"]),
                latest_analysis_at=values["latest_analysis_at"],
                latest_therapy_at=values["latest_therapy_at"],
                databases=tuple(sorted(values["databases"], key=str)),
            )
            for gene, values in sorted(aggregate.items())
        )

    def snapshots(self, gene: str, kind: str) -> tuple[ExecutionSnapshot, ...]:
        """Role: Load analysis or therapy snapshots for one gene."""
        if kind not in KNOWN_KINDS:
            raise ValueError(f"Unknown execution kind: {kind}")
        normalized_gene = gene.strip().upper()
        snapshots: list[ExecutionSnapshot] = []
        for path in self.available_databases():
            if not _has_table(path, "execution_cache"):
                continue
            with closing(_connect_read_only(path)) as connection:
                rows = connection.execute(
                    """
                    SELECT id, kind, gene, parameter_hash, parameters_json,
                           result_json, created_at, last_accessed_at, cache_hit_count
                    FROM execution_cache
                    WHERE gene = ? AND kind = ?
                    ORDER BY created_at DESC, id DESC
                    """,
                    (normalized_gene, kind),
                ).fetchall()
            for row in rows:
                try:
                    parameters = json.loads(str(row["parameters_json"]))
                    result = json.loads(str(row["result_json"]))
                except (TypeError, json.JSONDecodeError) as exc:
                    raise CatalogError(
                        f"Invalid JSON in {path}, execution {row['id']}: {exc}"
                    ) from exc
                if not isinstance(parameters, dict) or not isinstance(result, dict):
                    raise CatalogError(
                        f"Execution {row['id']} in {path} does not contain JSON objects"
                    )
                snapshots.append(
                    ExecutionSnapshot(
                        database_path=path,
                        execution_id=int(row["id"]),
                        kind=str(row["kind"]),
                        gene=str(row["gene"]),
                        parameter_hash=str(row["parameter_hash"]),
                        parameters=parameters,
                        result=result,
                        created_at=str(row["created_at"]),
                        last_accessed_at=str(row["last_accessed_at"]),
                        cache_hit_count=int(row["cache_hit_count"]),
                    )
                )
        snapshots.sort(
            key=lambda snapshot: (
                _parse_iso(snapshot.created_at),
                snapshot.execution_id,
                str(snapshot.database_path),
            ),
            reverse=True,
        )
        return tuple(snapshots)

    def latest(self, gene: str, kind: str) -> ExecutionSnapshot | None:
        """Role: Return the latest matching snapshot for a gene and result kind."""
        snapshots = self.snapshots(gene, kind)
        return snapshots[0] if snapshots else None

    def string_networks(self, gene: str) -> tuple[StringNetworkSnapshot, ...]:
        """Role: Load stored STRING-network snapshots for one gene."""
        normalized_gene = gene.strip().upper()
        snapshots: list[StringNetworkSnapshot] = []
        for path in self.available_databases():
            if not _has_table(path, "string_network_cache"):
                continue
            with closing(_connect_read_only(path)) as connection:
                rows = connection.execute(
                    """
                    SELECT gene, species, string_version, required_score,
                           network_type, added_nodes, result_json, created_at,
                           last_accessed_at, cache_hit_count
                    FROM string_network_cache
                    WHERE gene = ?
                    ORDER BY created_at DESC
                    """,
                    (normalized_gene,),
                ).fetchall()
            for row in rows:
                try:
                    result = json.loads(str(row["result_json"]))
                except (TypeError, json.JSONDecodeError) as exc:
                    raise CatalogError(
                        f"Invalid STRING JSON in {path} for {normalized_gene}: {exc}"
                    ) from exc
                if not isinstance(result, dict):
                    raise CatalogError(
                        f"STRING cache in {path} does not contain a JSON object"
                    )
                snapshots.append(
                    StringNetworkSnapshot(
                        database_path=path,
                        gene=str(row["gene"]),
                        species=int(row["species"]),
                        string_version=str(row["string_version"]),
                        required_score=int(row["required_score"]),
                        network_type=str(row["network_type"]),
                        added_nodes=int(row["added_nodes"]),
                        result=result,
                        created_at=str(row["created_at"]),
                        last_accessed_at=str(row["last_accessed_at"]),
                        cache_hit_count=int(row["cache_hit_count"]),
                    )
                )
        snapshots.sort(
            key=lambda snapshot: (
                _parse_iso(snapshot.created_at),
                snapshot.required_score,
                snapshot.added_nodes,
                str(snapshot.database_path),
            ),
            reverse=True,
        )
        return tuple(snapshots)

    def latest_string_network(self, gene: str) -> StringNetworkSnapshot | None:
        """Role: Return the latest stored STRING network for one gene."""
        snapshots = self.string_networks(gene)
        return snapshots[0] if snapshots else None


def _connect_read_only(path: Path) -> sqlite3.Connection:
    """Role: Compute connect read only for the catalog workflow."""
    uri = f"{path.resolve().as_uri()}?mode=ro"
    connection = sqlite3.connect(uri, uri=True, timeout=2.0)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA query_only = ON")
    return connection


def _has_table(path: Path, table: str) -> bool:
    """Role: Compute has table for the catalog workflow."""
    if not path.is_file():
        return False
    try:
        with closing(_connect_read_only(path)) as connection:
            row = connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
                (table,),
            ).fetchone()
        return row is not None
    except sqlite3.Error:
        return False


def _is_supported_database(path: Path) -> bool:
    """Role: Compute is supported database for the catalog workflow."""
    return _has_table(path, "execution_cache") or _has_table(
        path, "string_network_cache"
    )


def _parse_iso(value: str) -> datetime:
    """Role: Parse and validate iso from its source representation."""
    normalized = value.replace("Z", "+00:00")
    parsed = datetime.fromisoformat(normalized)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _newest_iso(current: str | None, candidate: str) -> str:
    """Role: Compute newest iso for the catalog workflow."""
    if current is None or _parse_iso(candidate) > _parse_iso(current):
        return candidate
    return current
