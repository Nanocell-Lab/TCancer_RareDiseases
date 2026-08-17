"""Generate validated mutant-receptor inputs for AlphaFold 3 and Boltz-2."""

from __future__ import annotations

import hashlib
import json
import os
import re
import sqlite3
from collections.abc import Iterable
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, Protocol

from .models import AnalysisResult
from .normalization import normalize_protein_change
from .protein_annotations import ProteinAnnotation, UniProtProteinClient


_CANONICAL_AA = "ACDEFGHIKLMNPQRSTVWY"
_MISSENSE = re.compile(rf"^([{_CANONICAL_AA}])(\d+)([{_CANONICAL_AA}])$")


class StructureInputError(ValueError):
    """Raised when a trustworthy mutant sequence cannot be generated."""


class ProteinSequenceClient(Protocol):
    def fetch(self, gene: str) -> ProteinAnnotation | None:
        """Role: Retrieve the reference protein sequence for model generation."""
        ...


@dataclass(frozen=True)
class StructureInputArtifact:
    platform: str
    gene: str
    protein_change: str
    accession: str | None
    sequence_length: int
    sequence_sha256: str
    path: str


@dataclass(frozen=True)
class SkippedStructureVariant:
    protein_change: str
    reason: str


@dataclass(frozen=True)
class StructureInputResult:
    gene: str
    accession: str | None
    sequence_length: int
    artifacts: tuple[StructureInputArtifact, ...]
    skipped: tuple[SkippedStructureVariant, ...]
    manifest_path: str

    def to_dict(self) -> dict[str, object]:
        """Role: Convert the domain object into a JSON-safe dictionary."""
        return asdict(self)


class ProteinSequenceCache:
    """Small cache for reviewed canonical UniProt sequences."""

    def __init__(self, path: str | Path) -> None:
        """Role: Compute init for the structure_inputs workflow."""
        if str(path) == ":memory:":
            raise ValueError("ProteinSequenceCache requires a persistent path")
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
        with self._connect() as connection:
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS protein_sequence_cache (
                    gene TEXT PRIMARY KEY,
                    accession TEXT,
                    sequence TEXT NOT NULL,
                    source_url TEXT,
                    retrieved_at TEXT NOT NULL,
                    last_accessed_at TEXT NOT NULL,
                    cache_hit_count INTEGER NOT NULL DEFAULT 0
                )
                """
            )
        self._restrict_permissions()

    def load(self, gene: str) -> tuple[str, str | None] | None:
        """Role: Load a cached resource when it is available and valid."""
        normalized = gene.strip().upper()
        with self._connect() as connection:
            row = connection.execute(
                """
                SELECT sequence, accession
                FROM protein_sequence_cache
                WHERE gene = ?
                """,
                (normalized,),
            ).fetchone()
            if row is None:
                return None
            connection.execute(
                """
                UPDATE protein_sequence_cache
                SET last_accessed_at = ?, cache_hit_count = cache_hit_count + 1
                WHERE gene = ?
                """,
                (_utc_now(), normalized),
            )
        return _validated_sequence(str(row["sequence"])), (
            str(row["accession"]) if row["accession"] else None
        )

    def store(self, annotation: ProteinAnnotation) -> None:
        """Role: Persist a normalized resource in the local cache."""
        if not annotation.sequence:
            raise StructureInputError("An empty sequence cannot be cached")
        sequence = _validated_sequence(annotation.sequence)
        now = _utc_now()
        with self._connect() as connection:
            connection.execute(
                """
                INSERT INTO protein_sequence_cache (
                    gene, accession, sequence, source_url, retrieved_at,
                    last_accessed_at, cache_hit_count
                ) VALUES (?, ?, ?, ?, ?, ?, 0)
                ON CONFLICT(gene) DO UPDATE SET
                    accession = excluded.accession,
                    sequence = excluded.sequence,
                    source_url = excluded.source_url,
                    retrieved_at = excluded.retrieved_at,
                    last_accessed_at = excluded.last_accessed_at,
                    cache_hit_count = 0
                """,
                (
                    annotation.gene.strip().upper(),
                    annotation.accession,
                    sequence,
                    annotation.source_url,
                    now,
                    now,
                ),
            )
        self._restrict_permissions()

    @contextmanager
    def _connect(self) -> Iterator[sqlite3.Connection]:
        """Role: Compute connect for the structure_inputs workflow."""
        connection = sqlite3.connect(self.path)
        connection.row_factory = sqlite3.Row
        try:
            with connection:
                yield connection
        finally:
            connection.close()

    def _restrict_permissions(self) -> None:
        """Role: Compute restrict permissions for the structure_inputs workflow."""
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


def generate_variant_structure_inputs(
    result: AnalysisResult,
    output_dir: str | Path,
    *,
    variants: Iterable[str] | None = None,
    protein_sequence: str | None = None,
    protein_accession: str | None = None,
    protein_client: ProteinSequenceClient | None = None,
    sequence_database: str | Path | None = None,
    offline: bool = False,
    model_seeds: Iterable[int] = (1,),
    max_variants: int = 20,
) -> StructureInputResult:
    """Write one AF3 JSON and one Boltz-2 YAML for each validated missense.

    When ``variants`` is omitted, shared rare-disease/oncology mutations are
    selected. The mutation is applied only when its reference residue agrees
    with the supplied or reviewed canonical UniProt sequence.
    """

    if max_variants < 1:
        raise ValueError("max_variants must be positive")
    seeds = _model_seeds(model_seeds)
    sequence, accession = _resolve_sequence(
        result.gene,
        protein_sequence=protein_sequence,
        protein_accession=protein_accession,
        protein_client=protein_client,
        sequence_database=sequence_database,
        offline=offline,
    )
    selected = _selected_variants(result, variants)[:max_variants]
    if not selected:
        raise StructureInputError(
            "No shared or explicit variants are available for structure generation"
        )

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    artifacts: list[StructureInputArtifact] = []
    skipped: list[SkippedStructureVariant] = []
    for change in selected:
        parsed = _MISSENSE.fullmatch(change)
        if parsed is None:
            skipped.append(
                SkippedStructureVariant(
                    protein_change=change,
                    reason="only simple missense substitutions are supported",
                )
            )
            continue
        reference, position_text, alternate = parsed.groups()
        position = int(position_text)
        if position > len(sequence):
            skipped.append(
                SkippedStructureVariant(
                    protein_change=change,
                    reason=(
                        f"position {position} is outside the canonical sequence "
                        f"of {len(sequence)} residues"
                    ),
                )
            )
            continue
        observed = sequence[position - 1]
        if observed != reference:
            skipped.append(
                SkippedStructureVariant(
                    protein_change=change,
                    reason=(
                        f"reference residue mismatch: expected "
                        f"{reference}{position}, but UniProt contains {observed}{position}"
                    ),
                )
            )
            continue
        mutant = sequence[: position - 1] + alternate + sequence[position:]
        stem = f"{result.gene.lower()}_{change.lower()}"
        description = (
            f"{result.gene} canonical receptor with p.{change}"
            + (f"; UniProt {accession}" if accession else "")
        )
        af3_path = output / f"{stem}_alphafold3.json"
        af3_path.write_text(
            json.dumps(
                _alphafold3_payload(
                    stem, mutant, description=description, model_seeds=seeds
                ),
                ensure_ascii=False,
                indent=2,
                allow_nan=False,
            )
            + "\n",
            encoding="utf-8",
        )
        boltz_path = output / f"{stem}_boltz2.yaml"
        boltz_path.write_text(
            _boltz2_yaml(
                mutant,
                gene=result.gene,
                protein_change=change,
                accession=accession,
            ),
            encoding="utf-8",
        )
        sequence_hash = hashlib.sha256(mutant.encode("ascii")).hexdigest()
        artifacts.extend(
            (
                StructureInputArtifact(
                    platform="alphafold3",
                    gene=result.gene,
                    protein_change=change,
                    accession=accession,
                    sequence_length=len(mutant),
                    sequence_sha256=sequence_hash,
                    path=str(af3_path),
                ),
                StructureInputArtifact(
                    platform="boltz2",
                    gene=result.gene,
                    protein_change=change,
                    accession=accession,
                    sequence_length=len(mutant),
                    sequence_sha256=sequence_hash,
                    path=str(boltz_path),
                ),
            )
        )

    manifest_path = output / f"{result.gene.lower()}_structure_inputs.json"
    generated = StructureInputResult(
        gene=result.gene,
        accession=accession,
        sequence_length=len(sequence),
        artifacts=tuple(artifacts),
        skipped=tuple(skipped),
        manifest_path=str(manifest_path),
    )
    manifest_path.write_text(
        json.dumps(
            generated.to_dict(),
            ensure_ascii=False,
            indent=2,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )
    return generated


def read_protein_sequence(path: str | Path) -> str:
    """Read a plain one-letter sequence or a single-record FASTA file."""

    source = Path(path)
    if not source.is_file():
        raise StructureInputError(f"Sequence file does not exist: {source}")
    lines = source.read_text(encoding="utf-8").splitlines()
    if sum(line.startswith(">") for line in lines) > 1:
        raise StructureInputError("The FASTA file must contain a single sequence")
    sequence = "".join(
        line.strip() for line in lines if line.strip() and not line.startswith(">")
    )
    return _validated_sequence(sequence)


def _resolve_sequence(
    gene: str,
    *,
    protein_sequence: str | None,
    protein_accession: str | None,
    protein_client: ProteinSequenceClient | None,
    sequence_database: str | Path | None,
    offline: bool,
) -> tuple[str, str | None]:
    """Role: Resolve sequence from cached, supplied, or inferred data."""
    if protein_sequence is not None:
        return _validated_sequence(protein_sequence), protein_accession
    cache = (
        ProteinSequenceCache(sequence_database)
        if sequence_database is not None
        else None
    )
    cached = cache.load(gene) if cache is not None else None
    if cached is not None:
        return cached
    if offline:
        raise StructureInputError(
            f"No UniProt sequence is cached for {gene}; "
            "provide --protein-sequence-file"
        )
    annotation = (protein_client or UniProtProteinClient()).fetch(gene)
    if annotation is None or not annotation.sequence:
        raise StructureInputError(
            f"UniProt did not return a reviewed canonical sequence for {gene}"
        )
    sequence = _validated_sequence(annotation.sequence)
    if annotation.length != len(sequence):
        raise StructureInputError(
            "The UniProt-reported length does not match its sequence"
        )
    if cache is not None:
        cache.store(annotation)
    return sequence, annotation.accession


def _validated_sequence(value: str) -> str:
    """Role: Compute validated sequence for the structure_inputs workflow."""
    sequence = re.sub(r"\s+", "", value).upper()
    if not sequence:
        raise StructureInputError("The protein sequence is empty")
    invalid = sorted(set(sequence) - set(_CANONICAL_AA))
    if invalid:
        raise StructureInputError(
            "The sequence contains non-canonical residues: " + ", ".join(invalid)
        )
    return sequence


def _selected_variants(
    result: AnalysisResult, variants: Iterable[str] | None
) -> tuple[str, ...]:
    """Role: Compute selected variants for the structure_inputs workflow."""
    raw = (
        tuple(variants)
        if variants is not None
        else tuple(item.protein_change for item in result.shared_mutations)
    )
    normalized: list[str] = []
    for value in raw:
        change = normalize_protein_change(value)
        if change and change not in normalized:
            normalized.append(change)
    return tuple(normalized)


def _model_seeds(values: Iterable[int]) -> tuple[int, ...]:
    """Role: Compute model seeds for the structure_inputs workflow."""
    seeds = tuple(values)
    if not seeds:
        raise ValueError("AlphaFold 3 requires at least one seed")
    if any(isinstance(seed, bool) or seed < 0 for seed in seeds):
        raise ValueError("Seeds must be zero or positive integers")
    return seeds


def _alphafold3_payload(
    name: str,
    sequence: str,
    *,
    description: str,
    model_seeds: tuple[int, ...],
) -> dict[str, object]:
    """Role: Compute alphafold3 payload for the structure_inputs workflow."""
    return {
        "name": name,
        "modelSeeds": list(model_seeds),
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": sequence,
                    "description": description,
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 4,
    }


def _boltz2_yaml(
    sequence: str,
    *,
    gene: str,
    protein_change: str,
    accession: str | None,
) -> str:
    """Role: Compute boltz2 yaml for the structure_inputs workflow."""
    comments = [
        f"# gene: {gene}",
        f"# variant: p.{protein_change}",
        f"# canonical_accession: {accession or 'not_provided'}",
        "# Run with: boltz predict INPUT.yaml --use_msa_server",
    ]
    return "\n".join(
        (
            *comments,
            "version: 1",
            "sequences:",
            "  - protein:",
            "      id: A",
            f"      sequence: {json.dumps(sequence)}",
            "",
        )
    )


def _utc_now() -> str:
    """Role: Return the current UTC timestamp in canonical form."""
    return datetime.now(timezone.utc).isoformat()
