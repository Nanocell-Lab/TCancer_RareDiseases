"""AlphaMissense annotations through Ensembl VEP with variant-level SQLite cache."""

from __future__ import annotations

import json
import math
import os
import sqlite3
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Protocol

from .evidence import EvidenceContext, gdc_identity, spdi_identity
from .http import ApiError, JsonHttpClient
from .models import EvidenceRecord
from .normalization import normalize_protein_change


ENSEMBL_VEP_REGION = "https://rest.ensembl.org/vep/homo_sapiens/region"
ALPHAMISSENSE_VERSION = "AlphaMissense-2023-GRCh38"
ALPHAMISSENSE_LICENSE = "CC BY 4.0"
COMPUTATIONAL_PREDICTION = "computational_prediction"


class AlphaMissenseHttp(Protocol):
    def post_json(self, url: str, payload: dict[str, Any]) -> Any:
        """Role: Submit variants to the configured annotation endpoint."""
        ...


@dataclass(frozen=True)
class AlphaMissensePrediction:
    gene: str
    variant_key: str
    assembly: str
    chromosome: str
    position: int
    reference: str
    alternate: str
    transcript: str
    protein_change: str
    score: float
    prediction_class: str
    canonical: bool
    mane_select: str | None
    retrieved_at: str


@dataclass(frozen=True)
class _Candidate:
    variant_key: str
    assembly: str
    chromosome: str
    position: int
    reference: str
    alternate: str
    protein_changes: tuple[str, ...]
    origins: tuple[str, ...]


@dataclass(frozen=True)
class _CacheEntry:
    predictions: tuple[AlphaMissensePrediction, ...]


class AlphaMissenseCache:
    """Persistent per-gene/per-allele cache stored beside execution snapshots."""

    def __init__(self, path: str | Path) -> None:
        """Role: Compute init for the alphamissense workflow."""
        if str(path) == ":memory:":
            raise ValueError("AlphaMissenseCache requires a persistent SQLite path")
        resolved = Path(path).expanduser().resolve()
        resolved.parent.mkdir(parents=True, exist_ok=True)
        if os.name == "posix" and not resolved.exists():
            try:
                descriptor = os.open(
                    resolved, os.O_CREAT | os.O_EXCL | os.O_RDWR, 0o600
                )
            except FileExistsError:
                pass
            else:
                os.close(descriptor)
        self.path = resolved
        self._initialize()

    def load(
        self,
        gene: str,
        variant_key: str,
        *,
        dataset_version: str = ALPHAMISSENSE_VERSION,
    ) -> _CacheEntry | None:
        """Role: Load a cached resource when it is available and valid."""
        with self._connect() as connection:
            row = connection.execute(
                """
                SELECT predictions_json
                FROM alphamissense_variant_cache
                WHERE gene = ? AND variant_key = ? AND dataset_version = ?
                """,
                (gene.strip().upper(), variant_key, dataset_version),
            ).fetchone()
            if row is None:
                return None
            connection.execute(
                """
                UPDATE alphamissense_variant_cache
                SET last_accessed_at = ?, cache_hit_count = cache_hit_count + 1
                WHERE gene = ? AND variant_key = ? AND dataset_version = ?
                """,
                (
                    _utc_now(),
                    gene.strip().upper(),
                    variant_key,
                    dataset_version,
                ),
            )
        payload = json.loads(str(row["predictions_json"]))
        if not isinstance(payload, list):
            raise ValueError("Invalid AlphaMissense SQLite record")
        return _CacheEntry(
            predictions=tuple(AlphaMissensePrediction(**item) for item in payload)
        )

    def store(
        self,
        gene: str,
        variant_key: str,
        predictions: tuple[AlphaMissensePrediction, ...],
        *,
        dataset_version: str = ALPHAMISSENSE_VERSION,
    ) -> None:
        """Role: Persist a normalized resource in the local cache."""
        now = _utc_now()
        payload = json.dumps(
            [asdict(item) for item in predictions],
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        with self._connect() as connection:
            connection.execute(
                """
                INSERT INTO alphamissense_variant_cache (
                    gene, variant_key, dataset_version, predictions_json,
                    created_at, last_accessed_at, cache_hit_count
                ) VALUES (?, ?, ?, ?, ?, ?, 0)
                ON CONFLICT(gene, variant_key, dataset_version) DO UPDATE SET
                    predictions_json = excluded.predictions_json,
                    created_at = excluded.created_at,
                    last_accessed_at = excluded.last_accessed_at,
                    cache_hit_count = 0
                """,
                (
                    gene.strip().upper(),
                    variant_key,
                    dataset_version,
                    payload,
                    now,
                    now,
                ),
            )
        self._restrict_permissions()

    @contextmanager
    def _connect(self) -> Iterator[sqlite3.Connection]:
        """Role: Compute connect for the alphamissense workflow."""
        connection = sqlite3.connect(self.path)
        connection.row_factory = sqlite3.Row
        try:
            with connection:
                yield connection
        finally:
            connection.close()

    def _initialize(self) -> None:
        """Role: Compute initialize for the alphamissense workflow."""
        with self._connect() as connection:
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS alphamissense_variant_cache (
                    gene TEXT NOT NULL,
                    variant_key TEXT NOT NULL,
                    dataset_version TEXT NOT NULL,
                    predictions_json TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    last_accessed_at TEXT NOT NULL,
                    cache_hit_count INTEGER NOT NULL DEFAULT 0,
                    PRIMARY KEY (gene, variant_key, dataset_version)
                )
                """
            )
        self._restrict_permissions()

    def _restrict_permissions(self) -> None:
        """Role: Compute restrict permissions for the alphamissense workflow."""
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


class AlphaMissenseEvidenceProvider:
    """Annotate exact GRCh38 SNVs and reuse prior variant results from SQLite."""

    name = "alphamissense"
    evidence_layer = COMPUTATIONAL_PREDICTION
    url: str | None = ENSEMBL_VEP_REGION
    source_version: str | None = ALPHAMISSENSE_VERSION

    def __init__(
        self,
        *,
        database: str | Path | None = None,
        http: AlphaMissenseHttp | None = None,
        batch_size: int = 200,
    ) -> None:
        """Role: Compute init for the alphamissense workflow."""
        if batch_size < 1 or batch_size > 200:
            raise ValueError("batch_size must be between 1 and 200")
        self.batch_size = batch_size
        self.cache = AlphaMissenseCache(database) if database is not None else None
        self.http = http or JsonHttpClient(
            user_agent="rare-oncology/2.1", min_interval=0.1
        )
        self.last_status = "ok"
        self.last_message: str | None = None

    def cache_identity(self) -> dict[str, Any]:
        """Role: Return the stable provider configuration used for cache keys."""
        return {
            "url": ENSEMBL_VEP_REGION,
            "dataset_version": ALPHAMISSENSE_VERSION,
            "batch_size": self.batch_size,
            "matching": "exact_grch38_snv_gene_and_protein_change",
        }

    def collect(self, context: EvidenceContext) -> list[EvidenceRecord]:
        """Role: Collect source-specific records for the current analysis context."""
        candidates, skipped = _candidates(context)
        cached: dict[str, tuple[AlphaMissensePrediction, ...]] = {}
        missing: list[_Candidate] = []
        cache_hits = 0
        for candidate in candidates:
            entry = (
                self.cache.load(context.gene, candidate.variant_key)
                if self.cache is not None
                else None
            )
            if entry is None:
                missing.append(candidate)
            else:
                cached[candidate.variant_key] = entry.predictions
                cache_hits += 1

        web_queries = 0
        for start in range(0, len(missing), self.batch_size):
            batch = missing[start : start + self.batch_size]
            queried = self._query_batch(context.gene, batch)
            web_queries += len(batch)
            for candidate in batch:
                predictions = queried.get(candidate.variant_key, ())
                cached[candidate.variant_key] = predictions
                if self.cache is not None:
                    self.cache.store(
                        context.gene, candidate.variant_key, predictions
                    )

        records: list[EvidenceRecord] = []
        unavailable = 0
        protein_mismatch = 0
        for candidate in candidates:
            predictions = cached.get(candidate.variant_key, ())
            matching = tuple(
                prediction
                for prediction in predictions
                if prediction.protein_change in candidate.protein_changes
            )
            if not predictions:
                unavailable += 1
            elif not matching:
                protein_mismatch += 1
            records.extend(
                _as_evidence(prediction, candidate) for prediction in matching
            )

        messages = [
            f"{cache_hits} variants from SQLite",
            f"{web_queries} variants queried in Ensembl VEP",
            f"{len(records)} predictions available",
        ]
        if unavailable:
            messages.append(f"{unavailable} without AlphaMissense data")
        if protein_mismatch:
            messages.append(
                f"{protein_mismatch} no exact protein match"
            )
        if skipped:
            messages.append(
                f"{skipped} changes cannot be queried (not an unambiguous GRCh38 SNV)"
            )
        self.last_status = "ok"
        self.last_message = "; ".join(messages)
        return records

    def _query_batch(
        self, gene: str, candidates: list[_Candidate]
    ) -> dict[str, tuple[AlphaMissensePrediction, ...]]:
        """Role: Compute query batch for the alphamissense workflow."""
        id_to_candidate = {
            f"am{index}": candidate for index, candidate in enumerate(candidates)
        }
        variants = [
            (
                f"{candidate.chromosome} {candidate.position} {identifier} "
                f"{candidate.reference} {candidate.alternate} . . ."
            )
            for identifier, candidate in id_to_candidate.items()
        ]
        response = self.http.post_json(
            ENSEMBL_VEP_REGION
            + "?AlphaMissense=1&canonical=1&mane=1&symbol=1",
            {"variants": variants},
        )
        if not isinstance(response, list):
            raise ApiError("Ensembl VEP returned a non-tabular JSON payload")
        collected: dict[str, list[AlphaMissensePrediction]] = {
            candidate.variant_key: [] for candidate in candidates
        }
        retrieved_at = _utc_now()
        for item in response:
            if not isinstance(item, Mapping):
                continue
            candidate = id_to_candidate.get(str(item.get("id", "")))
            if candidate is None:
                continue
            consequences = item.get("transcript_consequences") or []
            if not isinstance(consequences, list):
                continue
            for consequence in consequences:
                prediction = _parse_prediction(
                    gene, candidate, consequence, retrieved_at
                )
                if prediction is not None:
                    collected[candidate.variant_key].append(prediction)
        return {
            key: tuple(
                sorted(
                    {item.transcript: item for item in predictions}.values(),
                    key=lambda item: item.transcript,
                )
            )
            for key, predictions in collected.items()
        }


def _candidates(context: EvidenceContext) -> tuple[tuple[_Candidate, ...], int]:
    """Role: Compute candidates for the alphamissense workflow."""
    rows: dict[str, dict[str, Any]] = {}
    skipped = 0
    source_rows = [
        (spdi_identity(variant.canonical_spdi), variant.protein_change, "clinvar")
        for variant in context.rare_variants
    ]
    source_rows.extend(
        (gdc_identity(item.genomic_dna_change), item.protein_change, "gdc")
        for item in context.cancer_observations
    )
    for identity, protein_change, origin in source_rows:
        if not _exact_snv(identity):
            skipped += 1
            continue
        key = str(identity["variant_key"])
        row = rows.setdefault(
            key,
            {
                **identity,
                "protein_changes": set(),
                "origins": set(),
            },
        )
        normalized = normalize_protein_change(protein_change)
        if normalized:
            row["protein_changes"].add(normalized)
        row["origins"].add(origin)
    candidates = tuple(
        _Candidate(
            variant_key=key,
            assembly=str(row["assembly"]),
            chromosome=str(row["chromosome"]),
            position=int(row["position"]),
            reference=str(row["reference"]),
            alternate=str(row["alternate"]),
            protein_changes=tuple(sorted(row["protein_changes"])),
            origins=tuple(sorted(row["origins"])),
        )
        for key, row in sorted(rows.items())
    )
    return candidates, skipped


def _exact_snv(identity: Mapping[str, Any]) -> bool:
    """Role: Compute exact snv for the alphamissense workflow."""
    reference = str(identity.get("reference", "")).upper()
    alternate = str(identity.get("alternate", "")).upper()
    return (
        identity.get("assembly") == "GRCh38"
        and bool(identity.get("variant_key"))
        and len(reference) == 1
        and len(alternate) == 1
        and reference in "ACGT"
        and alternate in "ACGT"
    )


def _parse_prediction(
    gene: str,
    candidate: _Candidate,
    consequence: Any,
    retrieved_at: str,
) -> AlphaMissensePrediction | None:
    """Role: Parse and validate prediction from its source representation."""
    if not isinstance(consequence, Mapping):
        return None
    if str(consequence.get("gene_symbol", "")).upper() != gene.upper():
        return None
    terms = consequence.get("consequence_terms") or []
    if "missense_variant" not in terms:
        return None
    alphamissense = consequence.get("alphamissense")
    if not isinstance(alphamissense, Mapping):
        return None
    try:
        score = float(alphamissense["am_pathogenicity"])
        protein_position = int(consequence["protein_start"])
    except (KeyError, TypeError, ValueError):
        return None
    if not math.isfinite(score) or not 0 <= score <= 1:
        return None
    amino_acids = str(consequence.get("amino_acids", ""))
    reference_aa, separator, alternate_aa = amino_acids.partition("/")
    protein_change = normalize_protein_change(
        f"{reference_aa}{protein_position}{alternate_aa}"
        if separator
        else None
    )
    transcript = str(consequence.get("transcript_id", "")).strip()
    prediction_class = str(alphamissense.get("am_class", "")).strip()
    if not protein_change or not transcript or not prediction_class:
        return None
    mane_select = consequence.get("mane_select")
    return AlphaMissensePrediction(
        gene=gene.upper(),
        variant_key=candidate.variant_key,
        assembly=candidate.assembly,
        chromosome=candidate.chromosome,
        position=candidate.position,
        reference=candidate.reference,
        alternate=candidate.alternate,
        transcript=transcript,
        protein_change=protein_change,
        score=score,
        prediction_class=prediction_class,
        canonical=bool(consequence.get("canonical")),
        mane_select=str(mane_select) if mane_select else None,
        retrieved_at=retrieved_at,
    )


def _as_evidence(
    prediction: AlphaMissensePrediction, candidate: _Candidate
) -> EvidenceRecord:
    """Role: Compute as evidence for the alphamissense workflow."""
    return EvidenceRecord(
        source="alphamissense",
        source_record_id=(
            f"{prediction.variant_key}|{prediction.transcript}"
        ),
        evidence_layer=COMPUTATIONAL_PREDICTION,
        gene=prediction.gene,
        retrieved_at=prediction.retrieved_at,
        source_version=ALPHAMISSENSE_VERSION,
        license=ALPHAMISSENSE_LICENSE,
        variant_key=prediction.variant_key,
        assembly=prediction.assembly,
        chromosome=prediction.chromosome,
        position=prediction.position,
        reference=prediction.reference,
        alternate=prediction.alternate,
        transcript=prediction.transcript,
        protein_change=prediction.protein_change,
        evidence_level="computational_prediction",
        prediction_score=prediction.score,
        prediction_class=prediction.prediction_class,
        details=tuple(
            sorted(
                {
                    "model": "AlphaMissense",
                    "match_method": "exact_genomic_snv_gene_and_protein_change",
                    "origins": ",".join(candidate.origins),
                    "canonical_transcript": str(prediction.canonical).lower(),
                    "mane_select": prediction.mane_select or "",
                    "clinical_use": "not_validated_or_approved",
                }.items()
            )
        ),
    )


def _utc_now() -> str:
    """Role: Return the current UTC timestamp in canonical form."""
    return datetime.now(timezone.utc).isoformat()
