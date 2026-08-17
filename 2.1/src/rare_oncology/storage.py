"""SQLite execution history and cache-first orchestration."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import sqlite3
import urllib.parse
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable, Generic, Mapping, TypeVar

from .analysis import find_shared_mutations
from .models import (
    AnalysisResult,
    EvidenceRecord,
    MutationPlotRecord,
    SharedMutation,
    SourceStatus,
    TherapyAssociation,
    TherapyResult,
)
from .therapy import TherapyProvider, find_therapy_associations


SCHEMA_VERSION = 1
ANALYSIS_KIND = "mutation_analysis"
THERAPY_KIND = "therapy_analysis"
_KINDS = {ANALYSIS_KIND, THERAPY_KIND}

T = TypeVar("T", AnalysisResult, TherapyResult)


class CacheMissError(RuntimeError):
    """Raised when offline mode cannot satisfy a request from SQLite."""


@dataclass(frozen=True)
class CachedRun(Generic[T]):
    result: T
    cache_hit: bool
    execution_id: int
    stored_at: str
    database_path: str


@dataclass(frozen=True)
class ExecutionSummary:
    execution_id: int
    kind: str
    gene: str
    created_at: str
    last_accessed_at: str
    cache_hit_count: int
    parameter_hash: str


@dataclass(frozen=True)
class _StoredPayload:
    execution_id: int
    created_at: str
    result: AnalysisResult | TherapyResult


class ExecutionDatabase:
    """Append-only execution snapshots with indexed cache lookup."""

    def __init__(
        self,
        path: str | Path,
        *,
        clock: Callable[[], datetime] | None = None,
    ) -> None:
        """Role: Compute init for the storage workflow."""
        self.path = str(path)
        if self.path != ":memory:":
            database_path = Path(self.path).expanduser().resolve()
            database_path.parent.mkdir(parents=True, exist_ok=True)
            if os.name == "posix" and not database_path.exists():
                try:
                    descriptor = os.open(
                        database_path,
                        os.O_CREAT | os.O_EXCL | os.O_RDWR,
                        0o600,
                    )
                except FileExistsError:
                    pass
                else:
                    os.close(descriptor)
            self.path = str(database_path)
        self._clock = clock or (lambda: datetime.now(timezone.utc))
        self._connection = sqlite3.connect(self.path)
        self._connection.row_factory = sqlite3.Row
        self._connection.execute("PRAGMA foreign_keys = ON")
        if self.path != ":memory:":
            self._connection.execute("PRAGMA journal_mode = WAL")
        self._initialize()
        self._restrict_permissions()

    def close(self) -> None:
        """Role: Close the owned database connection."""
        self._connection.close()

    def __enter__(self) -> "ExecutionDatabase":
        """Role: Compute enter for the storage workflow."""
        return self

    def __exit__(self, *_exc: object) -> None:
        """Role: Compute exit for the storage workflow."""
        self.close()

    def save(
        self,
        kind: str,
        gene: str,
        parameters: Mapping[str, Any],
        result: AnalysisResult | TherapyResult,
    ) -> int:
        """Role: Persist a validated execution snapshot in SQLite."""
        _validate_kind(kind)
        normalized_gene = gene.strip().upper()
        expected_type = AnalysisResult if kind == ANALYSIS_KIND else TherapyResult
        if not isinstance(result, expected_type):
            raise TypeError(
                f"{kind} requires {expected_type.__name__}, not "
                f"{type(result).__name__}"
            )
        if result.gene.strip().upper() != normalized_gene:
            raise ValueError(
                f"The execution gene ({normalized_gene}) does not match "
                f"the result ({result.gene})"
            )
        parameters_json, parameter_hash = _parameter_identity(parameters)
        result_json = json.dumps(
            result.to_dict(),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        now = _as_utc(self._clock()).isoformat()
        with self._connection:
            cursor = self._connection.execute(
                """
                INSERT INTO execution_cache (
                    kind, gene, parameter_hash, parameters_json, result_json,
                    created_at, last_accessed_at, cache_hit_count, schema_version
                ) VALUES (?, ?, ?, ?, ?, ?, ?, 0, ?)
                """,
                (
                    kind,
                    normalized_gene,
                    parameter_hash,
                    parameters_json,
                    result_json,
                    now,
                    now,
                    SCHEMA_VERSION,
                ),
            )
        if cursor.lastrowid is None:
            raise RuntimeError("SQLite did not return an execution identifier")
        self._restrict_permissions()
        return int(cursor.lastrowid)

    def load_latest(
        self,
        kind: str,
        gene: str,
        parameters: Mapping[str, Any],
        *,
        max_age_hours: float | None = 24.0,
    ) -> _StoredPayload | None:
        """Role: Load the newest compatible execution snapshot from SQLite."""
        _validate_kind(kind)
        if max_age_hours is not None and (
            not math.isfinite(max_age_hours) or max_age_hours < 0
        ):
            raise ValueError("max_age_hours must be finite, zero, or positive")
        _, parameter_hash = _parameter_identity(parameters)
        row = self._connection.execute(
            """
            SELECT id, created_at, result_json
            FROM execution_cache
            WHERE kind = ? AND gene = ? AND parameter_hash = ?
              AND schema_version = ?
            ORDER BY id DESC
            LIMIT 1
            """,
            (kind, gene.strip().upper(), parameter_hash, SCHEMA_VERSION),
        ).fetchone()
        if row is None:
            return None
        created_at = _parse_datetime(str(row["created_at"]))
        if max_age_hours is not None:
            age = _as_utc(self._clock()) - created_at
            if age > timedelta(hours=max_age_hours):
                return None
        payload = json.loads(str(row["result_json"]))
        result = (
            analysis_result_from_dict(payload)
            if kind == ANALYSIS_KIND
            else therapy_result_from_dict(payload)
        )
        accessed_at = _as_utc(self._clock()).isoformat()
        with self._connection:
            self._connection.execute(
                """
                UPDATE execution_cache
                SET last_accessed_at = ?, cache_hit_count = cache_hit_count + 1
                WHERE id = ?
                """,
                (accessed_at, int(row["id"])),
            )
        self._restrict_permissions()
        return _StoredPayload(
            execution_id=int(row["id"]),
            created_at=created_at.isoformat(),
            result=result,
        )

    def list_executions(
        self,
        *,
        gene: str | None = None,
        kind: str | None = None,
        limit: int = 100,
    ) -> tuple[ExecutionSummary, ...]:
        """Role: List stored execution metadata without reconstructing results."""
        if limit < 1:
            raise ValueError("limit must be positive")
        if kind is not None:
            _validate_kind(kind)
        clauses: list[str] = []
        values: list[Any] = []
        if gene:
            clauses.append("gene = ?")
            values.append(gene.strip().upper())
        if kind:
            clauses.append("kind = ?")
            values.append(kind)
        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        values.append(limit)
        rows = self._connection.execute(
            f"""
            SELECT id, kind, gene, created_at, last_accessed_at,
                   cache_hit_count, parameter_hash
            FROM execution_cache
            {where}
            ORDER BY id DESC
            LIMIT ?
            """,
            values,
        ).fetchall()
        return tuple(
            ExecutionSummary(
                execution_id=int(row["id"]),
                kind=str(row["kind"]),
                gene=str(row["gene"]),
                created_at=str(row["created_at"]),
                last_accessed_at=str(row["last_accessed_at"]),
                cache_hit_count=int(row["cache_hit_count"]),
                parameter_hash=str(row["parameter_hash"]),
            )
            for row in rows
        )

    def _initialize(self) -> None:
        """Role: Compute initialize for the storage workflow."""
        with self._connection:
            self._connection.execute(
                """
                CREATE TABLE IF NOT EXISTS execution_cache (
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
                )
                """
            )
            self._connection.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_execution_cache_lookup
                ON execution_cache (
                    kind, gene, parameter_hash, schema_version, id DESC
                )
                """
            )

    def _restrict_permissions(self) -> None:
        """Role: Compute restrict permissions for the storage workflow."""
        if os.name != "posix" or self.path == ":memory:":
            return
        for candidate in (self.path, f"{self.path}-wal", f"{self.path}-shm"):
            try:
                os.chmod(candidate, 0o600)
            except FileNotFoundError:
                continue


def find_shared_mutations_cached(
    gene: str,
    *,
    database: ExecutionDatabase | str | Path,
    cache_parameters: Mapping[str, Any] | None = None,
    max_age_hours: float | None = 24.0,
    refresh: bool = False,
    offline: bool = False,
    **analysis_kwargs: Any,
) -> CachedRun[AnalysisResult]:
    """Read a compatible molecular analysis from SQLite before using web APIs."""

    materialized_kwargs = dict(analysis_kwargs)
    for key in _ANALYSIS_ITERABLE_ARGUMENTS:
        value = materialized_kwargs.get(key)
        if value is not None and not isinstance(value, (str, bytes, bytearray)):
            materialized_kwargs[key] = tuple(value)
    identity_kwargs = dict(materialized_kwargs)
    evidence_providers = identity_kwargs.get("evidence_providers")
    if evidence_providers is not None:
        identity_kwargs["evidence_providers"] = _provider_identities(
            tuple(evidence_providers), "evidence"
        )
    parameters = _cache_parameters(cache_parameters, identity_kwargs)
    db, owned = _database(database)
    try:
        if not refresh:
            stored = db.load_latest(
                ANALYSIS_KIND,
                gene,
                parameters,
                max_age_hours=max_age_hours,
            )
            if stored is not None:
                if not isinstance(stored.result, AnalysisResult):
                    raise TypeError("Unexpected molecular snapshot type")
                return CachedRun(
                    result=stored.result,
                    cache_hit=True,
                    execution_id=stored.execution_id,
                    stored_at=stored.created_at,
                    database_path=db.path,
                )
        if offline:
            raise CacheMissError(
                f"No current molecular analysis is cached for {gene.upper()}"
            )
        result = find_shared_mutations(gene, **materialized_kwargs)
        execution_id = db.save(ANALYSIS_KIND, gene, parameters, result)
        stored_at = _as_utc(db._clock()).isoformat()
        return CachedRun(
            result=result,
            cache_hit=False,
            execution_id=execution_id,
            stored_at=stored_at,
            database_path=db.path,
        )
    finally:
        if owned:
            db.close()


def find_therapy_associations_cached(
    result: AnalysisResult,
    *,
    database: ExecutionDatabase | str | Path,
    providers: Iterable[TherapyProvider] | None = None,
    cache_parameters: Mapping[str, Any] | None = None,
    max_age_hours: float | None = 24.0,
    refresh: bool = False,
    offline: bool = False,
    strict: bool = False,
) -> CachedRun[TherapyResult]:
    """Read therapy layers from SQLite before complementing them with web APIs."""

    provider_tuple = tuple(providers) if providers is not None else None
    defaults: dict[str, Any] = {
        "analysis_signature": _analysis_signature(result),
        "providers": _provider_identities(provider_tuple, "therapy"),
        "strict": strict,
    }
    if cache_parameters is not None:
        defaults["request"] = cache_parameters
    parameters = _cache_parameters(None, defaults)
    db, owned = _database(database)
    try:
        if not refresh:
            stored = db.load_latest(
                THERAPY_KIND,
                result.gene,
                parameters,
                max_age_hours=max_age_hours,
            )
            if stored is not None:
                if not isinstance(stored.result, TherapyResult):
                    raise TypeError("Therapy snapshot has an unexpected type")
                return CachedRun(
                    result=stored.result,
                    cache_hit=True,
                    execution_id=stored.execution_id,
                    stored_at=stored.created_at,
                    database_path=db.path,
                )
        if offline:
            raise CacheMissError(
                f"No current therapy analysis is cached for "
                f"{result.gene}"
            )
        therapy_result = find_therapy_associations(
            result, providers=provider_tuple, strict=strict
        )
        execution_id = db.save(
            THERAPY_KIND, result.gene, parameters, therapy_result
        )
        stored_at = _as_utc(db._clock()).isoformat()
        return CachedRun(
            result=therapy_result,
            cache_hit=False,
            execution_id=execution_id,
            stored_at=stored_at,
            database_path=db.path,
        )
    finally:
        if owned:
            db.close()


def analysis_result_from_dict(payload: Mapping[str, Any]) -> AnalysisResult:
    """Role: Provide the analysis result from dict operation for the storage workflow."""
    data = dict(payload)
    data.pop("association_count", None)
    data["shared_mutations"] = tuple(
        SharedMutation(**_tuple_fields(item, _SHARED_TUPLES))
        for item in data.get("shared_mutations", [])
    )
    data["mutation_plot_data"] = tuple(
        MutationPlotRecord(**_tuple_fields(item, _PLOT_TUPLES))
        for item in data.get("mutation_plot_data", [])
    )
    data["evidence_records"] = tuple(
        EvidenceRecord(**_tuple_fields(item, _EVIDENCE_TUPLES, nested=("details",)))
        for item in data.get("evidence_records", [])
    )
    data["source_statuses"] = tuple(
        SourceStatus(**item) for item in data.get("source_statuses", [])
    )
    for field in (
        "gdc_projects",
        "excluded_conditions",
        "allowed_germline_classifications",
        "sources",
    ):
        if data.get(field) is not None:
            data[field] = tuple(data[field])
    return AnalysisResult(**data)


def therapy_result_from_dict(payload: Mapping[str, Any]) -> TherapyResult:
    """Role: Provide the therapy result from dict operation for the storage workflow."""
    data = dict(payload)
    data.pop("association_count", None)
    for layer in (
        "variant_therapy_evidence",
        "regulatory_therapy",
        "trial_and_target_discovery",
    ):
        data[layer] = tuple(
            TherapyAssociation(
                **_tuple_fields(item, _THERAPY_TUPLES, nested=("details",))
            )
            for item in data.get(layer, [])
        )
    data["protein_changes"] = tuple(data.get("protein_changes", []))
    data["source_statuses"] = tuple(
        SourceStatus(**item) for item in data.get("source_statuses", [])
    )
    return TherapyResult(**data)


_SHARED_TUPLES = {
    "diseases",
    "orphanet_codes",
    "clinvar_ids",
    "clinvar_accessions",
    "germline_classifications",
    "canonical_spdis",
    "tcga_projects",
    "genomic_confirmations",
    "ssm_ids",
}
_PLOT_TUPLES = {
    "diseases",
    "orphanet_codes",
    "germline_classifications",
    "project_patient_counts",
    "disease_variant_counts",
}
_EVIDENCE_TUPLES = {"details"}
_THERAPY_TUPLES = {
    "therapy_components",
    "disease_ids",
    "disease_names",
    "details",
}
_ANALYSIS_ITERABLE_ARGUMENTS = {
    "projects",
    "exclude_conditions",
    "allowed_germline_classifications",
    "evidence_providers",
}


def _tuple_fields(
    payload: Mapping[str, Any],
    fields: set[str],
    *,
    nested: tuple[str, ...] = (),
) -> dict[str, Any]:
    """Role: Compute tuple fields for the storage workflow."""
    data = dict(payload)
    for field in fields:
        if field in data and data[field] is not None:
            data[field] = tuple(data[field])
    for field in nested:
        if field in data:
            data[field] = tuple(tuple(item) for item in data[field])
    if "project_patient_counts" in data:
        data["project_patient_counts"] = tuple(
            tuple(item) for item in data["project_patient_counts"]
        )
    if "disease_variant_counts" in data:
        data["disease_variant_counts"] = tuple(
            tuple(item) for item in data["disease_variant_counts"]
        )
    return data


def _cache_parameters(
    explicit: Mapping[str, Any] | None,
    inferred: Mapping[str, Any],
) -> dict[str, Any]:
    """Role: Compute cache parameters for the storage workflow."""
    source: dict[str, Any] = {"arguments": inferred}
    if explicit is not None:
        source["request"] = explicit
    normalized = _normalize_cache_value(source)
    if not isinstance(normalized, dict):
        raise TypeError("Cache parameters must be a mapping")
    return normalized


def _normalize_cache_value(value: Any) -> Any:
    """Role: Normalize cache value into a stable internal representation."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        resolved = value.expanduser().resolve()
        descriptor: dict[str, Any] = {"path": str(resolved)}
        if resolved.exists():
            stat = resolved.stat()
            descriptor.update({"size": stat.st_size, "mtime_ns": stat.st_mtime_ns})
            if resolved.is_file():
                descriptor["sha256"] = _file_sha256(resolved)
        return descriptor
    if isinstance(value, Mapping):
        normalized_mapping: dict[str, Any] = {}
        for key, item in sorted(value.items(), key=lambda pair: str(pair[0])):
            name = str(key)
            folded = name.casefold()
            if _is_sensitive_name(folded) or folded == "email":
                normalized_mapping[name] = "REDACTED"
            elif (folded == "url" or folded.endswith("_url")) and isinstance(
                item, str
            ):
                normalized_mapping[name] = _redact_url(item)
            else:
                normalized_mapping[name] = _normalize_cache_value(item)
        return normalized_mapping
    if isinstance(value, (list, tuple, set, frozenset)):
        items = [_normalize_cache_value(item) for item in value]
        return (
            sorted(items, key=_canonical_json)
            if isinstance(value, (set, frozenset))
            else items
        )
    explicit_identity = getattr(value, "cache_identity", None)
    if callable(explicit_identity):
        return {
            "type": f"{type(value).__module__}.{type(value).__qualname__}",
            "configuration": _normalize_cache_value(explicit_identity()),
        }
    descriptor: dict[str, Any] = {
        "type": f"{type(value).__module__}.{type(value).__qualname__}"
    }
    for attribute in (
        "name",
        "layer",
        "evidence_layer",
        "max_variants",
        "max_studies",
        "max_drugs",
        "path",
        "jurisdiction",
        "include_details",
        "page_size",
        "source_version",
        "url",
    ):
        if hasattr(value, attribute):
            attribute_value = getattr(value, attribute)
            if attribute == "url" and isinstance(attribute_value, str):
                descriptor[attribute] = _redact_url(attribute_value)
            else:
                descriptor[attribute] = _normalize_cache_value(attribute_value)
    return descriptor


def _provider_identities(
    providers: tuple[Any, ...] | None,
    provider_kind: str,
) -> tuple[dict[str, Any], ...] | None:
    """Role: Resolve or describe the configured identities provider."""
    if providers is None:
        return None
    identities: list[dict[str, Any]] = []
    for provider in providers:
        identity_method = getattr(provider, "cache_identity", None)
        if not callable(identity_method):
            raise TypeError(
                f"The {provider_kind} provider {type(provider).__name__} must "
                "implement cache_identity() with stable configuration"
            )
        identity = identity_method()
        if not isinstance(identity, Mapping):
            raise TypeError("cache_identity() must return a mapping")
        identities.append(
            {
                "type": f"{type(provider).__module__}.{type(provider).__qualname__}",
                "configuration": dict(identity),
            }
        )
    return tuple(identities)


def _analysis_signature(result: AnalysisResult) -> dict[str, Any]:
    """Role: Compute analysis signature for the storage workflow."""
    return {
        "gene": result.gene,
        "shared_mutations": [
            {
                "protein_change": item.protein_change,
                "diseases": list(item.diseases),
                "orphanet_codes": list(item.orphanet_codes),
                "match_level": item.match_level,
            }
            for item in result.shared_mutations
        ],
    }


def _parameter_identity(parameters: Mapping[str, Any]) -> tuple[str, str]:
    """Role: Compute parameter identity for the storage workflow."""
    parameters_json = _canonical_json(_normalize_cache_value(parameters))
    return parameters_json, hashlib.sha256(parameters_json.encode("utf-8")).hexdigest()


def _canonical_json(value: Any) -> str:
    """Role: Compute canonical json for the storage workflow."""
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _database(
    database: ExecutionDatabase | str | Path,
) -> tuple[ExecutionDatabase, bool]:
    """Role: Compute database for the storage workflow."""
    if isinstance(database, ExecutionDatabase):
        return database, False
    return ExecutionDatabase(database), True


def _validate_kind(kind: str) -> None:
    """Role: Validate kind constraints before execution."""
    if kind not in _KINDS:
        raise ValueError(f"Unknown execution type: {kind}")


def _parse_datetime(value: str) -> datetime:
    """Role: Parse and validate datetime from its source representation."""
    return _as_utc(datetime.fromisoformat(value.replace("Z", "+00:00")))


def _as_utc(value: datetime) -> datetime:
    """Role: Compute as utc for the storage workflow."""
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _redact_url(value: str) -> str:
    """Role: Compute redact url for the storage workflow."""
    parsed = urllib.parse.urlsplit(value)
    query = urllib.parse.parse_qsl(parsed.query, keep_blank_values=True)
    redacted = [
        (key, "REDACTED" if _is_sensitive_name(key, url_context=True) else item)
        for key, item in query
    ]
    netloc = parsed.netloc.rsplit("@", 1)[-1]
    fragment = (
        "REDACTED"
        if _contains_sensitive_marker(parsed.fragment)
        else parsed.fragment
    )
    return urllib.parse.urlunsplit(
        (
            parsed.scheme,
            netloc,
            parsed.path,
            urllib.parse.urlencode(redacted),
            fragment,
        )
    )


def _file_sha256(path: Path) -> str:
    """Role: Compute or inspect sha256 file metadata."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _is_sensitive_name(value: str, *, url_context: bool = False) -> bool:
    """Role: Compute is sensitive name for the storage workflow."""
    normalized = re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_")
    parts = set(normalized.split("_"))
    return bool(
        parts.intersection(
            {"token", "secret", "password", "credential", "credentials"}
        )
        or (
            "signature" in parts
            and (
                url_context
                or bool(
                    parts.intersection(
                        {"aws", "amz", "auth", "request", "client"}
                    )
                )
            )
        )
        or normalized
        in {
            "key",
            "apikey",
            "api_key",
            "access_key",
            "secret_key",
            "private_key",
            "client_key",
            "consumer_key",
            "aws_access_key_id",
            "auth",
            "authorization",
            "bearer",
        }
    )


def _contains_sensitive_marker(value: str) -> bool:
    """Role: Compute contains sensitive marker for the storage workflow."""
    return any(
        marker in value.casefold()
        for marker in ("token", "secret", "password", "signature", "credential")
    )
