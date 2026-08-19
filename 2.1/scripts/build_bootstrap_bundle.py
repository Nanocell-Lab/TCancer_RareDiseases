#!/usr/bin/env python3
"""Build the tracked, compact dashboard bootstrap from full local snapshots."""

from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


GENES = ("EGFR", "TGFBR1", "ACVR1", "PDGFRB")
MAX_BOOTSTRAP_BYTES = 10 * 1024 * 1024
ANALYSIS_KIND = "mutation_analysis"
THERAPY_KIND = "therapy_analysis"
SCHEMA_VERSION = 1
ROOT = Path(__file__).resolve().parents[1]


EXECUTION_SCHEMA = """
CREATE TABLE execution_cache (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    kind TEXT NOT NULL,
    gene TEXT NOT NULL,
    parameter_hash TEXT NOT NULL,
    parameters_json TEXT NOT NULL,
    result_json TEXT NOT NULL,
    created_at TEXT NOT NULL,
    last_accessed_at TEXT NOT NULL,
    cache_hit_count INTEGER NOT NULL DEFAULT 0,
    schema_version INTEGER NOT NULL
);
CREATE INDEX idx_execution_cache_lookup
ON execution_cache (kind, gene, parameter_hash, schema_version, id DESC);
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
);
"""

METADATA_SCHEMA = """
CREATE TABLE bootstrap_metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
"""


def build_parser() -> argparse.ArgumentParser:
    """Define reproducible source and destination locations for the bundle."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, default=ROOT / "results")
    parser.add_argument(
        "--regulatory-root", type=Path, default=ROOT / "downloads" / "regulatory"
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "bootstrap")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Build the SQLite, copy normalized regulatory tables and write provenance."""

    args = build_parser().parse_args(argv)
    manifest = build_bundle(
        args.results_root.expanduser().resolve(),
        args.regulatory_root.expanduser().resolve(),
        args.output_dir.expanduser().resolve(),
    )
    print(
        f"Bootstrap: {manifest['genes']} -> {manifest['database']['path']} "
        f"({manifest['database']['size_bytes']} bytes)"
    )
    return 0


def build_bundle(
    results_root: Path,
    regulatory_root: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Create a bounded read-only seed while retaining dashboard projections."""

    databases = tuple(
        sorted(
            path
            for path in results_root.rglob("*.sqlite3")
            if path.is_file() and path.name != "global_family_index.sqlite3"
        )
    )
    if not databases:
        raise FileNotFoundError(f"No source SQLite databases found under {results_root}")

    analyses = _execution_rows(databases, ANALYSIS_KIND)
    therapies = _execution_rows(databases, THERAPY_KIND)
    networks = _string_rows(databases)
    selected_analyses = {
        gene: _newest(row for row in analyses if row["gene"] == gene)
        for gene in GENES
    }
    selected_therapies = {
        gene: _matching_therapy(selected_analyses[gene], therapies)
        for gene in GENES
    }
    selected_networks = {
        gene: _newest_or_none(row for row in networks if row["gene"] == gene)
        for gene in GENES
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    database_path = output_dir / "rare_oncology_demo.sqlite3"
    temporary_path = output_dir / ".rare_oncology_demo.sqlite3.tmp"
    temporary_path.unlink(missing_ok=True)
    _write_database(
        temporary_path,
        selected_analyses,
        selected_therapies,
        selected_networks,
        results_root,
    )
    if temporary_path.stat().st_size >= MAX_BOOTSTRAP_BYTES:
        temporary_path.unlink(missing_ok=True)
        raise ValueError("Bootstrap SQLite exceeds the 10 MiB repository limit")
    temporary_path.replace(database_path)

    regulatory_dir = output_dir / "regulatory"
    regulatory_dir.mkdir(parents=True, exist_ok=True)
    regulatory_files = []
    for filename in ("fda_biomarkers.csv", "ema_labels.csv"):
        source = regulatory_root / filename
        if not source.is_file():
            raise FileNotFoundError(f"Normalized regulatory table not found: {source}")
        destination = regulatory_dir / filename
        destination.write_bytes(source.read_bytes().replace(b"\r\n", b"\n"))
        regulatory_files.append(_file_descriptor(destination, output_dir))

    manifest: dict[str, Any] = {
        "format_version": 1,
        "purpose": "Compact read-only examples for a functional Git clone",
        "genes": list(GENES),
        "database": _file_descriptor(database_path, output_dir),
        "regulatory_tables": regulatory_files,
        "source_snapshots": [
            {
                "gene": gene,
                "analysis": _source_descriptor(selected_analyses[gene], results_root),
                "therapy": _source_descriptor(selected_therapies[gene], results_root),
                "string": (
                    _source_descriptor(network, results_root)
                    if (network := selected_networks[gene]) is not None
                    else None
                ),
            }
            for gene in GENES
        ],
        "compaction": {
            "analysis": (
                "All shared mutations and mutation_plot_data are retained. Detailed "
                "ClinVar/GDC evidence rows are omitted; AlphaMissense is deduplicated "
                "to displayed protein changes, while gnomAD and PanelApp rows remain."
            ),
            "therapy": "The newest molecular-signature-compatible result is retained.",
            "string": "The newest stored human network is retained.",
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def _execution_rows(
    databases: Iterable[Path], kind: str
) -> list[dict[str, Any]]:
    """Read supported immutable execution rows from every source database."""

    rows: list[dict[str, Any]] = []
    for path in databases:
        with sqlite3.connect(f"file:{path}?mode=ro", uri=True) as connection:
            connection.row_factory = sqlite3.Row
            if not _has_table(connection, "execution_cache"):
                continue
            for row in connection.execute(
                """
                SELECT id, kind, gene, parameter_hash, parameters_json, result_json,
                       created_at, last_accessed_at, schema_version
                FROM execution_cache
                WHERE kind = ? AND schema_version = ?
                """,
                (kind, SCHEMA_VERSION),
            ):
                rows.append(
                    {
                        **dict(row),
                        "gene": str(row["gene"]).upper(),
                        "source_path": path,
                        "parameters": json.loads(str(row["parameters_json"])),
                        "result": json.loads(str(row["result_json"])),
                    }
                )
    return rows


def _string_rows(databases: Iterable[Path]) -> list[dict[str, Any]]:
    """Read stored STRING networks without contacting STRING."""

    rows: list[dict[str, Any]] = []
    for path in databases:
        with sqlite3.connect(f"file:{path}?mode=ro", uri=True) as connection:
            connection.row_factory = sqlite3.Row
            if not _has_table(connection, "string_network_cache"):
                continue
            for row in connection.execute(
                """
                SELECT gene, species, string_version, required_score, network_type,
                       added_nodes, result_json, created_at, last_accessed_at
                FROM string_network_cache
                """
            ):
                rows.append(
                    {
                        **dict(row),
                        "gene": str(row["gene"]).upper(),
                        "source_path": path,
                        "result": json.loads(str(row["result_json"])),
                    }
                )
    return rows


def _newest(rows: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Select the newest deterministic candidate or fail visibly."""

    candidates = list(rows)
    if not candidates:
        raise ValueError("A required bootstrap snapshot is missing")
    return max(
        candidates,
        key=lambda row: (
            _timestamp(str(row["created_at"])),
            int(row.get("id", 0)),
            str(row["source_path"]),
        ),
    )


def _newest_or_none(rows: Iterable[dict[str, Any]]) -> dict[str, Any] | None:
    """Select an optional newest snapshot without weakening required records."""

    candidates = list(rows)
    return _newest(candidates) if candidates else None


def _matching_therapy(
    analysis: Mapping[str, Any], therapies: Iterable[dict[str, Any]]
) -> dict[str, Any]:
    """Select therapy evidence created from exactly the same molecular signature."""

    signature = _analysis_signature(analysis["result"])
    return _newest(
        row
        for row in therapies
        if row["gene"] == analysis["gene"]
        and _therapy_signature(row["parameters"]) == signature
    )


def _analysis_signature(result: Mapping[str, Any]) -> dict[str, Any]:
    """Reproduce the stable signature stored in therapy cache parameters."""

    return {
        "gene": str(result.get("gene") or "").upper(),
        "shared_mutations": [
            {
                "protein_change": record.get("protein_change"),
                "diseases": list(record.get("diseases") or []),
                "orphanet_codes": list(record.get("orphanet_codes") or []),
                "match_level": record.get("match_level"),
            }
            for record in result.get("shared_mutations", [])
            if isinstance(record, dict)
        ],
    }


def _therapy_signature(parameters: Mapping[str, Any]) -> Mapping[str, Any] | None:
    """Extract the molecular signature without retaining machine-local paths."""

    arguments = parameters.get("arguments")
    if isinstance(arguments, dict) and isinstance(
        arguments.get("analysis_signature"), dict
    ):
        return arguments["analysis_signature"]
    return None


def _compact_analysis(result: Mapping[str, Any]) -> dict[str, Any]:
    """Keep every plotted aggregate while bounding detailed evidence volume."""

    compact = dict(result)
    changes = {
        str(record.get("protein_change"))
        for record in compact.get("mutation_plot_data", [])
        if isinstance(record, dict) and record.get("protein_change")
    }
    retained: list[Mapping[str, Any]] = []
    seen_alpha: set[tuple[Any, ...]] = set()
    for record in compact.get("evidence_records", []):
        if not isinstance(record, dict):
            continue
        source = str(record.get("source") or "").casefold()
        if source == "alphamissense":
            if str(record.get("protein_change")) not in changes:
                continue
            key = (
                record.get("protein_change"),
                record.get("prediction_score"),
                record.get("prediction_class"),
            )
            if key in seen_alpha:
                continue
            seen_alpha.add(key)
            retained.append(record)
        elif source in {"gnomad", "panelapp"}:
            retained.append(record)
    original_count = int(compact.get("evidence_record_count") or 0)
    compact["evidence_records"] = retained
    compact["evidence_record_count"] = len(retained)
    compact["caveat"] = (
        str(compact.get("caveat") or "").rstrip()
        + " Bootstrap snapshot: detailed ClinVar/GDC evidence rows were compacted "
        + f"from {original_count} total records; shared mutations, plot projections "
        + "and source provenance remain available."
    )
    return compact


def _compact_therapy(result: Mapping[str, Any]) -> dict[str, Any]:
    """Remove machine-local table paths while preserving therapy assertions."""

    compact = dict(result)
    statuses = []
    for value in compact.get("source_statuses", []):
        if not isinstance(value, dict):
            continue
        status = dict(value)
        source = str(status.get("source") or "").casefold()
        if source in {"fda", "ema"}:
            filename = (
                "fda_biomarkers.csv" if source == "fda" else "ema_labels.csv"
            )
            status["url"] = f"bootstrap/regulatory/{filename}"
        statuses.append(status)
    compact["source_statuses"] = statuses
    compact["caveat"] = (
        str(compact.get("caveat") or "").rstrip()
        + " Bootstrap snapshot: regulatory table locations are repository-relative."
    )
    return compact


def _write_database(
    path: Path,
    analyses: Mapping[str, Mapping[str, Any]],
    therapies: Mapping[str, Mapping[str, Any]],
    networks: Mapping[str, Mapping[str, Any] | None],
    results_root: Path,
) -> None:
    """Write a deterministic SQLite with only dashboard-supported public tables."""

    with sqlite3.connect(path) as connection:
        connection.executescript(EXECUTION_SCHEMA + STRING_SCHEMA + METADATA_SCHEMA)
        for gene in GENES:
            analysis = analyses[gene]
            analysis_parameters = {
                "bootstrap": True,
                "source_parameter_hash": analysis["parameter_hash"],
            }
            therapy = therapies[gene]
            therapy_parameters = {
                "arguments": {
                    "analysis_signature": _analysis_signature(analysis["result"])
                },
                "bootstrap": True,
                "source_parameter_hash": therapy["parameter_hash"],
            }
            for kind, source, parameters, result in (
                (
                    ANALYSIS_KIND,
                    analysis,
                    analysis_parameters,
                    _compact_analysis(analysis["result"]),
                ),
                (
                    THERAPY_KIND,
                    therapy,
                    therapy_parameters,
                    _compact_therapy(therapy["result"]),
                ),
            ):
                parameters_json = _canonical_json(parameters)
                connection.execute(
                    """
                    INSERT INTO execution_cache (
                        kind, gene, parameter_hash, parameters_json, result_json,
                        created_at, last_accessed_at, cache_hit_count, schema_version
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, 0, ?)
                    """,
                    (
                        kind,
                        gene,
                        hashlib.sha256(parameters_json.encode()).hexdigest(),
                        parameters_json,
                        _canonical_json(result),
                        source["created_at"],
                        source["last_accessed_at"],
                        SCHEMA_VERSION,
                    ),
                )
            network = networks[gene]
            if network is not None:
                connection.execute(
                    """
                    INSERT INTO string_network_cache (
                        gene, species, string_version, required_score, network_type,
                        added_nodes, result_json, created_at, last_accessed_at,
                        cache_hit_count
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, 0)
                    """,
                    (
                        gene,
                        network["species"],
                        network["string_version"],
                        network["required_score"],
                        network["network_type"],
                        network["added_nodes"],
                        _canonical_json(network["result"]),
                        network["created_at"],
                        network["last_accessed_at"],
                    ),
                )
        metadata = {
            "format_version": "1",
            "genes": ",".join(GENES),
            "source_root_label": results_root.name,
            "research_use_only": "true",
        }
        connection.executemany(
            "INSERT INTO bootstrap_metadata (key, value) VALUES (?, ?)",
            sorted(metadata.items()),
        )
        connection.commit()
        connection.execute("VACUUM")


def _source_descriptor(row: Mapping[str, Any], root: Path) -> dict[str, Any]:
    """Record portable provenance without exposing the builder's absolute path."""

    source = Path(row["source_path"])
    try:
        label = str(source.relative_to(root))
    except ValueError:
        label = source.name
    return {
        "database": label,
        "execution_id": row.get("id"),
        "created_at": row["created_at"],
    }


def _file_descriptor(path: Path, root: Path) -> dict[str, Any]:
    """Describe one tracked artifact for clone-time integrity checks."""

    return {
        "path": str(path.relative_to(root)),
        "size_bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def _has_table(connection: sqlite3.Connection, table: str) -> bool:
    """Return whether a source SQLite exposes the requested public table."""

    return (
        connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)
        ).fetchone()
        is not None
    )


def _timestamp(value: str) -> datetime:
    """Parse stored ISO timestamps into a comparable UTC value."""

    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    return parsed.replace(tzinfo=parsed.tzinfo or timezone.utc).astimezone(timezone.utc)


def _canonical_json(value: object) -> str:
    """Serialize cache payloads deterministically and reject non-finite numbers."""

    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


if __name__ == "__main__":
    raise SystemExit(main())
