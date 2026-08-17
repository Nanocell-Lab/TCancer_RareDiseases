"""Human STRING functional-association networks with SQLite reuse."""

from __future__ import annotations

import csv
import json
import math
import os
import re
import sqlite3
import urllib.parse
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Protocol

from .http import ApiError, JsonHttpClient


STRING_VERSION = "12.0"
STRING_BASE_URL = "https://version-12-0.string-db.org"
STRING_API_URL = f"{STRING_BASE_URL}/api/json"
STRING_SPECIES_HUMAN = 9606
STRING_CALLER_IDENTITY = "rare-oncology-2.1"
STRING_CAVEAT = (
    "STRING represents functional or physical protein associations according to "
    "network_type; an edge alone does not establish physical binding, causality, "
    "or therapeutic relevance."
)
_GENE_PATTERN = re.compile(r"^[A-Z][A-Z0-9.-]{0,39}$")
_NETWORK_TYPES = {"functional", "physical"}


class StringNetworkError(RuntimeError):
    """Raised when a STRING network cannot be resolved or reused."""


class StringHttp(Protocol):
    def get_json(self, url: str, params: dict[str, Any]) -> Any:
        """Role: Retrieve a JSON response from a STRING endpoint."""
        ...


@dataclass(frozen=True)
class StringNode:
    string_id: str
    preferred_name: str
    is_query: bool
    annotation: str | None = None


@dataclass(frozen=True)
class StringEdge:
    string_id_a: str
    string_id_b: str
    preferred_name_a: str
    preferred_name_b: str
    score: float
    neighborhood_score: float
    fusion_score: float
    phylogenetic_score: float
    coexpression_score: float
    experimental_score: float
    database_score: float
    textmining_score: float


@dataclass(frozen=True)
class StringNetwork:
    gene: str
    species: int
    string_version: str
    mapped_string_id: str
    mapped_preferred_name: str
    mapping_annotation: str | None
    required_score: int
    network_type: str
    added_nodes: int
    retrieved_at: str
    nodes: tuple[StringNode, ...]
    edges: tuple[StringEdge, ...]
    source_url: str
    caveat: str = STRING_CAVEAT

    def to_dict(self) -> dict[str, Any]:
        """Role: Convert the domain object into a JSON-safe dictionary."""
        return asdict(self)


@dataclass(frozen=True)
class StringNetworkRun:
    network: StringNetwork
    cache_hit: bool
    database_path: str | None


class StringNetworkCache:
    """Exact-parameter cache stored beside the 2.1 execution history."""

    def __init__(self, path: str | Path) -> None:
        """Role: Compute init for the string_network workflow."""
        if str(path) == ":memory:":
            raise ValueError("StringNetworkCache requires a persistent SQLite path")
        self.path = Path(path).expanduser().resolve()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if os.name == "posix" and not self.path.exists():
            try:
                descriptor = os.open(
                    self.path, os.O_CREAT | os.O_EXCL | os.O_RDWR, 0o600
                )
            except FileExistsError:
                pass
            else:
                os.close(descriptor)
        self._initialize()

    def load(
        self,
        gene: str,
        *,
        required_score: int,
        network_type: str,
        added_nodes: int,
    ) -> StringNetwork | None:
        """Role: Load a cached resource when it is available and valid."""
        with self._connect() as connection:
            row = connection.execute(
                """
                SELECT result_json
                FROM string_network_cache
                WHERE gene = ? AND species = ? AND string_version = ?
                  AND required_score = ? AND network_type = ? AND added_nodes = ?
                """,
                (
                    _normalize_gene(gene),
                    STRING_SPECIES_HUMAN,
                    STRING_VERSION,
                    required_score,
                    network_type,
                    added_nodes,
                ),
            ).fetchone()
            if row is None:
                return None
            connection.execute(
                """
                UPDATE string_network_cache
                SET last_accessed_at = ?, cache_hit_count = cache_hit_count + 1
                WHERE gene = ? AND species = ? AND string_version = ?
                  AND required_score = ? AND network_type = ? AND added_nodes = ?
                """,
                (
                    _utc_now(),
                    _normalize_gene(gene),
                    STRING_SPECIES_HUMAN,
                    STRING_VERSION,
                    required_score,
                    network_type,
                    added_nodes,
                ),
            )
        payload = json.loads(str(row["result_json"]))
        if not isinstance(payload, dict):
            raise StringNetworkError("Invalid STRING SQLite record")
        return string_network_from_dict(payload)

    def store(self, network: StringNetwork) -> None:
        """Role: Persist a normalized resource in the local cache."""
        result_json = json.dumps(
            network.to_dict(),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        now = _utc_now()
        with self._connect() as connection:
            connection.execute(
                """
                INSERT INTO string_network_cache (
                    gene, species, string_version, required_score, network_type,
                    added_nodes, result_json, created_at, last_accessed_at,
                    cache_hit_count
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, 0)
                ON CONFLICT(
                    gene, species, string_version, required_score, network_type,
                    added_nodes
                ) DO UPDATE SET
                    result_json = excluded.result_json,
                    created_at = excluded.created_at,
                    last_accessed_at = excluded.last_accessed_at,
                    cache_hit_count = 0
                """,
                (
                    network.gene,
                    network.species,
                    network.string_version,
                    network.required_score,
                    network.network_type,
                    network.added_nodes,
                    result_json,
                    now,
                    now,
                ),
            )
        self._restrict_permissions()

    @contextmanager
    def _connect(self) -> Iterator[sqlite3.Connection]:
        """Role: Compute connect for the string_network workflow."""
        connection = sqlite3.connect(self.path)
        connection.row_factory = sqlite3.Row
        try:
            with connection:
                yield connection
        finally:
            connection.close()

    def _initialize(self) -> None:
        """Role: Compute initialize for the string_network workflow."""
        with self._connect() as connection:
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS string_network_cache (
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
                        gene, species, string_version, required_score,
                        network_type, added_nodes
                    )
                )
                """
            )
        self._restrict_permissions()

    def _restrict_permissions(self) -> None:
        """Role: Compute restrict permissions for the string_network workflow."""
        if os.name != "posix":
            return
        for candidate in (
            self.path,
            Path(f"{self.path}-wal"),
            Path(f"{self.path}-shm"),
        ):
            try:
                os.chmod(candidate, 0o600)
            except FileNotFoundError:
                continue


def get_string_network(
    gene: str,
    *,
    database: str | Path | None = None,
    required_score: int = 400,
    network_type: str = "functional",
    added_nodes: int = 20,
    refresh: bool = False,
    offline: bool = False,
    http: StringHttp | None = None,
) -> StringNetworkRun:
    """Load an exact human STRING network from SQLite before querying STRING."""

    normalized_gene = _normalize_gene(gene)
    _validate_parameters(required_score, network_type, added_nodes)
    if refresh and offline:
        raise ValueError("refresh and offline cannot be combined")
    cache = StringNetworkCache(database) if database is not None else None
    if cache is not None and not refresh:
        cached = cache.load(
            normalized_gene,
            required_score=required_score,
            network_type=network_type,
            added_nodes=added_nodes,
        )
        if cached is not None:
            return StringNetworkRun(cached, True, str(cache.path))
    if offline:
        raise StringNetworkError(
            f"No compatible STRING network is cached in SQLite for {normalized_gene}"
        )
    client = http or JsonHttpClient(
        user_agent="rare-oncology/2.1",
        min_interval=1.0,
    )
    network = _fetch_string_network(
        normalized_gene,
        client,
        required_score=required_score,
        network_type=network_type,
        added_nodes=added_nodes,
    )
    if cache is not None:
        cache.store(network)
    return StringNetworkRun(network, False, str(cache.path) if cache else None)


def _fetch_string_network(
    gene: str,
    http: StringHttp,
    *,
    required_score: int,
    network_type: str,
    added_nodes: int,
) -> StringNetwork:
    """Role: Compute fetch string network for the string_network workflow."""
    mapping_url = f"{STRING_API_URL}/get_string_ids"
    common = {
        "species": STRING_SPECIES_HUMAN,
        "caller_identity": STRING_CALLER_IDENTITY,
    }
    try:
        mapping_payload = http.get_json(
            mapping_url,
            {"identifiers": gene, "echo_query": 1, **common},
        )
    except ApiError as exc:
        raise StringNetworkError(f"STRING could not map {gene}: {exc}") from exc
    mapping = _choose_mapping(gene, mapping_payload)
    string_id = _required_text(mapping, "stringId")
    if not string_id.startswith(f"{STRING_SPECIES_HUMAN}."):
        raise StringNetworkError(
            f"STRING resolved {gene} outside Homo sapiens: {string_id}"
        )
    preferred_name = _required_text(mapping, "preferredName")
    annotation = _optional_text(mapping.get("annotation"))
    network_url = f"{STRING_API_URL}/network"
    try:
        network_payload = http.get_json(
            network_url,
            {
                "identifiers": string_id,
                "required_score": required_score,
                "network_type": network_type,
                "add_nodes": added_nodes,
                **common,
            },
        )
    except ApiError as exc:
        raise StringNetworkError(
            f"STRING could not retrieve the network for {gene}: {exc}"
        ) from exc
    edges = _parse_edges(network_payload)
    nodes: dict[str, StringNode] = {
        string_id: StringNode(
            string_id=string_id,
            preferred_name=preferred_name,
            is_query=True,
            annotation=annotation,
        )
    }
    for edge in edges:
        for node_id, name in (
            (edge.string_id_a, edge.preferred_name_a),
            (edge.string_id_b, edge.preferred_name_b),
        ):
            existing = nodes.get(node_id)
            is_query = node_id == string_id
            nodes[node_id] = StringNode(
                string_id=node_id,
                preferred_name=(preferred_name if is_query else name),
                is_query=is_query,
                annotation=annotation if is_query else (existing.annotation if existing else None),
            )
    source_url = f"{STRING_BASE_URL}/cgi/network?" + urllib.parse.urlencode(
        {
            "identifiers": string_id,
            "species": STRING_SPECIES_HUMAN,
            "required_score": required_score,
            "network_type": network_type,
            "add_nodes": added_nodes,
        }
    )
    return StringNetwork(
        gene=gene,
        species=STRING_SPECIES_HUMAN,
        string_version=STRING_VERSION,
        mapped_string_id=string_id,
        mapped_preferred_name=preferred_name,
        mapping_annotation=annotation,
        required_score=required_score,
        network_type=network_type,
        added_nodes=added_nodes,
        retrieved_at=_utc_now(),
        nodes=tuple(
            sorted(nodes.values(), key=lambda node: (not node.is_query, node.preferred_name))
        ),
        edges=edges,
        source_url=source_url,
    )


def _choose_mapping(gene: str, payload: Any) -> dict[str, Any]:
    """Role: Compute choose mapping for the string_network workflow."""
    if not isinstance(payload, list) or not payload:
        raise StringNetworkError(f"STRING did not find a human identifier for {gene}")
    mappings = [item for item in payload if isinstance(item, dict)]
    if not mappings:
        raise StringNetworkError(f"Invalid STRING mapping response for {gene}")
    return next(
        (
            item
            for item in mappings
            if str(item.get("preferredName", "")).casefold() == gene.casefold()
        ),
        mappings[0],
    )


def _parse_edges(payload: Any) -> tuple[StringEdge, ...]:
    """Role: Parse and validate edges from its source representation."""
    if not isinstance(payload, list):
        raise StringNetworkError("Invalid STRING network response")
    deduplicated: dict[tuple[str, str], StringEdge] = {}
    for item in payload:
        if not isinstance(item, dict):
            continue
        id_a = _required_text(item, "stringId_A")
        id_b = _required_text(item, "stringId_B")
        if id_a == id_b:
            continue
        edge = StringEdge(
            string_id_a=id_a,
            string_id_b=id_b,
            preferred_name_a=_required_text(item, "preferredName_A"),
            preferred_name_b=_required_text(item, "preferredName_B"),
            score=_score(item, "score"),
            neighborhood_score=_score(item, "nscore"),
            fusion_score=_score(item, "fscore"),
            phylogenetic_score=_score(item, "pscore"),
            coexpression_score=_score(item, "ascore"),
            experimental_score=_score(item, "escore"),
            database_score=_score(item, "dscore"),
            textmining_score=_score(item, "tscore"),
        )
        key = (id_a, id_b) if id_a < id_b else (id_b, id_a)
        previous = deduplicated.get(key)
        if previous is None or edge.score > previous.score:
            deduplicated[key] = edge
    return tuple(
        sorted(
            deduplicated.values(),
            key=lambda edge: (-edge.score, edge.preferred_name_a, edge.preferred_name_b),
        )
    )


def string_network_from_dict(payload: dict[str, Any]) -> StringNetwork:
    """Role: Provide the string network from dict operation for the string_network workflow."""
    return StringNetwork(
        gene=str(payload["gene"]),
        species=int(payload["species"]),
        string_version=str(payload["string_version"]),
        mapped_string_id=str(payload["mapped_string_id"]),
        mapped_preferred_name=str(payload["mapped_preferred_name"]),
        mapping_annotation=_optional_text(payload.get("mapping_annotation")),
        required_score=int(payload["required_score"]),
        network_type=str(payload["network_type"]),
        added_nodes=int(payload["added_nodes"]),
        retrieved_at=str(payload["retrieved_at"]),
        nodes=tuple(StringNode(**node) for node in payload.get("nodes", [])),
        edges=tuple(StringEdge(**edge) for edge in payload.get("edges", [])),
        source_url=str(payload["source_url"]),
        caveat=str(payload.get("caveat") or STRING_CAVEAT),
    )


def write_string_network_json(network: StringNetwork, path: str | Path) -> Path:
    """Role: Serialize and write string network json output."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(network.to_dict(), ensure_ascii=False, indent=2, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )
    return output


def write_string_nodes_csv(network: StringNetwork, path: str | Path) -> Path:
    """Role: Serialize and write string nodes csv output."""
    rows = [asdict(node) for node in network.nodes]
    return _write_csv(path, rows, ("string_id", "preferred_name", "is_query", "annotation"))


def write_string_edges_csv(network: StringNetwork, path: str | Path) -> Path:
    """Role: Serialize and write string edges csv output."""
    rows = [asdict(edge) for edge in network.edges]
    fieldnames = tuple(StringEdge.__dataclass_fields__)
    return _write_csv(path, rows, fieldnames)


def _write_csv(
    path: str | Path,
    rows: list[dict[str, Any]],
    fieldnames: tuple[str, ...],
) -> Path:
    """Role: Compute write csv for the string_network workflow."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _safe_csv_value(value) for key, value in row.items()})
    return output


def _safe_csv_value(value: Any) -> Any:
    """Role: Sanitize csv value for safe serialization."""
    if not isinstance(value, str):
        return value
    trimmed = value.lstrip(" \t\r\n")
    if value.startswith(("\t", "\r", "\n")) or trimmed.startswith(("=", "+", "-", "@")):
        return "'" + value
    return value


def _score(item: dict[str, Any], key: str) -> float:
    """Role: Compute score for the string_network workflow."""
    raw = item.get(key, 0.0)
    try:
        value = float(raw)
    except (TypeError, ValueError) as exc:
        raise StringNetworkError(f"Invalid STRING score in {key}: {raw!r}") from exc
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise StringNetworkError(f"STRING score out of range in {key}: {raw!r}")
    return value


def _required_text(item: dict[str, Any], key: str) -> str:
    """Role: Compute required text for the string_network workflow."""
    value = _optional_text(item.get(key))
    if not value:
        raise StringNetworkError(f"Missing STRING field: {key}")
    return value


def _optional_text(value: Any) -> str | None:
    """Role: Convert optional text input into its typed value."""
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _normalize_gene(gene: str) -> str:
    """Role: Normalize gene into a stable internal representation."""
    normalized = gene.strip().upper()
    if not _GENE_PATTERN.fullmatch(normalized):
        raise ValueError(f"Invalid human gene symbol: {gene!r}")
    return normalized


def _validate_parameters(required_score: int, network_type: str, added_nodes: int) -> None:
    """Role: Validate parameters constraints before execution."""
    if not 0 <= required_score <= 1000:
        raise ValueError("required_score must be between 0 and 1000")
    if network_type not in _NETWORK_TYPES:
        raise ValueError("network_type must be 'functional' or 'physical'")
    if not 0 <= added_nodes <= 100:
        raise ValueError("added_nodes must be between 0 and 100")


def _utc_now() -> str:
    """Role: Return the current UTC timestamp in canonical form."""
    return datetime.now(timezone.utc).isoformat()
