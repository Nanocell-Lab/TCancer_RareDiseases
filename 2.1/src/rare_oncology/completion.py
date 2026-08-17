"""Backfill missing data-source layers for genes already stored in SQLite."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import sqlite3
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Protocol, Sequence

from .alphamissense import AlphaMissenseEvidenceProvider
from .evidence import (
    GRCH38_ACCESSIONS,
    EvidenceContext,
    EvidenceProvider,
    GnomadEvidenceProvider,
    PanelAppEvidenceProvider,
    collect_evidence,
)
from .models import (
    AnalysisResult,
    CancerObservation,
    EvidenceRecord,
    RareVariant,
    SourceStatus,
    TherapyAssociation,
    TherapyResult,
)
from .http import ApiError
from .protein_annotations import ProteinAnnotation, UniProtProteinClient
from .storage import (
    ANALYSIS_KIND,
    THERAPY_KIND,
    ExecutionDatabase,
    analysis_result_from_dict,
    therapy_result_from_dict,
)
from .string_network import (
    StringHttp,
    StringNetworkError,
    get_string_network,
)
from .structure_inputs import ProteinSequenceCache
from .therapy import (
    REGULATORY_THERAPY,
    TRIAL_AND_TARGET_DISCOVERY,
    VARIANT_THERAPY_EVIDENCE,
    CivicTherapyProvider,
    ClinicalTrialsTherapyProvider,
    OncoKBTherapyProvider,
    OpenTargetsTherapyProvider,
    RegulatoryTableTherapyProvider,
    TherapyContext,
    TherapyProvider,
    collect_therapy_associations,
)


DEFAULT_COMPLETION_SOURCES = (
    "string",
    "uniprot",
    "alphamissense",
    "gnomad",
    "panelapp",
    "civic",
    "clinicaltrials",
    "opentargets",
)
EVIDENCE_SOURCES = {"alphamissense", "gnomad", "panelapp"}
THERAPY_SOURCES = {"civic", "oncokb", "clinicaltrials", "opentargets"}
COMPLETE_STATUSES = {"ok", "partial", "unavailable"}
COMPLETION_SCHEMA_VERSION = 1


class ProteinClient(Protocol):
    def fetch(self, gene: str) -> ProteinAnnotation | None:
        """Role: Retrieve a normalized protein annotation for one gene."""
        ...


@dataclass(frozen=True)
class CompletionItem:
    database: str
    gene: str
    source: str
    action: str
    status: str
    record_count: int = 0
    message: str | None = None


@dataclass(frozen=True)
class CompletionReport:
    generated_at: str
    databases: tuple[str, ...]
    requested_sources: tuple[str, ...]
    gene_count: int
    query_count: int
    cache_hit_count: int
    completed_count: int
    missing_count: int
    error_count: int
    items: tuple[CompletionItem, ...]
    caveat: str = (
        "Each source retains its own layer and provenance; a STRING network, "
        "prediction, trial, or gene-drug association does not establish clinical actionability."
    )

    def to_dict(self) -> dict[str, Any]:
        """Role: Convert the domain object into a JSON-safe dictionary."""
        return asdict(self)


@dataclass(frozen=True)
class _AnalysisSnapshot:
    execution_id: int
    parameters: Mapping[str, Any]
    result: AnalysisResult


@dataclass(frozen=True)
class _TherapySnapshot:
    execution_id: int
    parameters: Mapping[str, Any]
    result: TherapyResult


def discover_completion_databases(
    results_root: str | Path, *, exclude: Iterable[str | Path] = ()
) -> tuple[Path, ...]:
    """Return SQLite files containing the 2.1 execution history table."""

    root = Path(results_root).expanduser().resolve()
    excluded = {Path(value).expanduser().resolve() for value in exclude}
    output: list[Path] = []
    for path in sorted(root.rglob("*.sqlite3"), key=str):
        resolved = path.resolve()
        if resolved in excluded or not resolved.is_file():
            continue
        try:
            with sqlite3.connect(f"{resolved.as_uri()}?mode=ro", uri=True) as connection:
                exists = connection.execute(
                    "SELECT 1 FROM sqlite_master WHERE type='table' AND name='execution_cache'"
                ).fetchone()
        except sqlite3.Error:
            continue
        if exists:
            output.append(resolved)
    return tuple(output)


def complete_missing_queries(
    databases: Iterable[str | Path],
    *,
    sources: Iterable[str] = DEFAULT_COMPLETION_SOURCES,
    oncokb_token: str | None = None,
    regulatory_files: Iterable[tuple[str, str | Path]] = (),
    max_therapy_records: int = 100,
    gnomad_max_variants: int = 25,
    string_score: int = 700,
    string_partners: int = 20,
    string_network_type: str = "functional",
    offline: bool = False,
    dry_run: bool = False,
    retry_errors: bool = False,
    strict: bool = False,
    max_genes: int | None = None,
    evidence_providers: Sequence[EvidenceProvider] | None = None,
    therapy_providers: Sequence[TherapyProvider] | None = None,
    string_http: StringHttp | None = None,
    protein_client: ProteinClient | None = None,
) -> CompletionReport:
    """Complete only missing source layers for every stored gene.

    Molecular ClinVar/GDC observations are never recomputed. Optional evidence is
    reconstructed from their stored exact identities and appended as a new
    compatible snapshot. Repeated executions first inspect the source-specific
    cache/status and therefore do not repeat completed queries.
    """

    regulatory_tuple = tuple(
        (source.strip().lower(), Path(path)) for source, path in regulatory_files
    )
    regulatory_names = {source for source, _path in regulatory_tuple}
    selected = tuple(
        dict.fromkeys(
            [value.strip().lower() for value in sources if value.strip()]
            + [source for source, _path in regulatory_tuple]
        )
    )
    if max_genes is not None and max_genes < 1:
        raise ValueError("max_genes must be positive")
    if max_therapy_records < 1:
        raise ValueError("max_therapy_records must be positive")
    unknown = set(selected) - (
        set(DEFAULT_COMPLETION_SOURCES) | {"oncokb"} | regulatory_names | {provider.name.lower() for provider in evidence_providers or ()} | {provider.name.lower() for provider in therapy_providers or ()}
    )
    if unknown:
        raise ValueError("Unknown completion sources: " + ", ".join(sorted(unknown)))
    if "oncokb" in selected and not oncokb_token and not any(
        provider.name.lower() == "oncokb" for provider in therapy_providers or ()
    ):
        raise ValueError("OncoKB requires oncokb_token/ONCOKB_TOKEN")

    paths = tuple(dict.fromkeys(Path(value).expanduser().resolve() for value in databases))
    items: list[CompletionItem] = []
    seen_genes = 0
    for database in paths:
        analyses, therapies = _read_snapshots(database)
        if not analyses:
            continue
        if not dry_run and not offline:
            _initialize_completion_status(database)
        for gene, snapshot in analyses.items():
            if max_genes is not None and seen_genes >= max_genes:
                break
            seen_genes += 1
            analysis = snapshot.result
            context = _evidence_context(analysis)
            analysis, evidence_items = _complete_evidence(
                database,
                snapshot.parameters,
                analysis,
                context,
                selected,
                offline=offline,
                dry_run=dry_run,
                retry_errors=retry_errors,
                strict=strict,
                gnomad_max_variants=gnomad_max_variants,
                injected=evidence_providers,
            )
            items.extend(evidence_items)
            therapy_rows = therapies.get(gene, ())
            therapy_items = _complete_therapy(
                database,
                analysis,
                therapy_rows,
                selected,
                oncokb_token=oncokb_token,
                regulatory_files=regulatory_tuple,
                max_records=max_therapy_records,
                offline=offline,
                dry_run=dry_run,
                retry_errors=retry_errors,
                strict=strict,
                injected=therapy_providers,
            )
            items.extend(therapy_items)
            if "string" in selected:
                items.append(
                    _complete_string(
                        database,
                        gene,
                        required_score=string_score,
                        added_nodes=string_partners,
                        network_type=string_network_type,
                        offline=offline,
                        dry_run=dry_run,
                        retry_errors=retry_errors,
                        http=string_http,
                    )
                )
            if "uniprot" in selected:
                items.append(
                    _complete_uniprot(
                        database,
                        gene,
                        offline=offline,
                        dry_run=dry_run,
                        retry_errors=retry_errors,
                        client=protein_client,
                    )
                )
        if max_genes is not None and seen_genes >= max_genes:
            break
    return CompletionReport(
        generated_at=_utc_now(),
        databases=tuple(str(path) for path in paths),
        requested_sources=selected,
        gene_count=len({(item.database, item.gene) for item in items}),
        query_count=sum(item.action == "queried" for item in items),
        cache_hit_count=sum(item.action == "cache_hit" for item in items),
        completed_count=sum(item.status in COMPLETE_STATUSES for item in items),
        missing_count=sum(item.status == "missing" for item in items),
        error_count=sum(item.status == "error" for item in items),
        items=tuple(items),
    )


def write_completion_report(
    report: CompletionReport, output_dir: str | Path
) -> tuple[Path, Path]:
    """Role: Serialize and write completion report output."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "completion_report.json"
    csv_path = output / "completion_report.csv"
    json_path.write_text(
        json.dumps(report.to_dict(), ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = tuple(CompletionItem.__dataclass_fields__)
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in report.items:
            writer.writerow({key: _safe_csv(value) for key, value in asdict(item).items()})
    return json_path, csv_path


def _read_snapshots(
    database: Path,
) -> tuple[dict[str, _AnalysisSnapshot], dict[str, tuple[_TherapySnapshot, ...]]]:
    """Role: Compute read snapshots for the completion workflow."""
    analyses: dict[str, _AnalysisSnapshot] = {}
    therapies: dict[str, list[_TherapySnapshot]] = {}
    try:
        with sqlite3.connect(f"{database.as_uri()}?mode=ro", uri=True) as connection:
            connection.row_factory = sqlite3.Row
            exists = connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='execution_cache'"
            ).fetchone()
            if not exists:
                return {}, {}
            rows = connection.execute(
                """
                SELECT id, kind, gene, parameters_json, result_json
                FROM execution_cache
                WHERE kind IN (?, ?)
                ORDER BY id DESC
                """,
                (ANALYSIS_KIND, THERAPY_KIND),
            ).fetchall()
    except sqlite3.Error:
        return {}, {}
    for row in rows:
        gene = str(row["gene"]).upper()
        payload = json.loads(str(row["result_json"]))
        if row["kind"] == ANALYSIS_KIND and gene not in analyses:
            analyses[gene] = _AnalysisSnapshot(
                int(row["id"]),
                json.loads(str(row["parameters_json"])),
                analysis_result_from_dict(payload),
            )
        elif row["kind"] == THERAPY_KIND:
            parameters = json.loads(str(row["parameters_json"]))
            therapies.setdefault(gene, []).append(
                _TherapySnapshot(
                    int(row["id"]), parameters, therapy_result_from_dict(payload)
                )
            )
    return analyses, {gene: tuple(rows) for gene, rows in therapies.items()}


def _complete_evidence(
    database: Path,
    base_parameters: Mapping[str, Any],
    analysis: AnalysisResult,
    context: EvidenceContext,
    selected: tuple[str, ...],
    *,
    offline: bool,
    dry_run: bool,
    retry_errors: bool,
    strict: bool,
    gnomad_max_variants: int,
    injected: Sequence[EvidenceProvider] | None,
) -> tuple[AnalysisResult, list[CompletionItem]]:
    """Role: Backfill missing evidence data for one stored gene."""
    requested = [source for source in selected if source in EVIDENCE_SOURCES]
    custom = {provider.name.lower(): provider for provider in injected or ()}
    requested.extend(source for source in selected if source in custom and source not in requested)
    if not requested:
        return analysis, []
    statuses = _latest_statuses(analysis.source_statuses)
    records = list(analysis.evidence_records)
    new_statuses: list[SourceStatus] = []
    items: list[CompletionItem] = []
    queried_sources: list[str] = []
    for source in requested:
        existing = statuses.get(source)
        if existing is not None and (
            _status_is_complete(existing)
            or (existing.status == "error" and not retry_errors)
        ):
            items.append(
                CompletionItem(
                    str(database), analysis.gene, source, "cache_hit",
                    existing.status, existing.record_count, existing.message,
                )
            )
            continue
        if dry_run or offline:
            items.append(
                CompletionItem(
                    str(database), analysis.gene, source,
                    "would_query" if dry_run else "offline_missing", "missing",
                    message="web query skipped",
                )
            )
            continue
        provider = custom.get(source) or _evidence_provider(
            source, database, gnomad_max_variants
        )
        identity = _provider_identity(provider)
        prior = _load_completion_status(database, analysis.gene, source, identity)
        if prior and (
            prior[0] in COMPLETE_STATUSES
            or (prior[0] == "error" and not retry_errors)
        ):
            items.append(
                CompletionItem(str(database), analysis.gene, source, "cache_hit", *prior)
            )
            continue
        collected, batch_statuses = collect_evidence(context, (provider,), strict=strict)
        status = batch_statuses[0]
        queried_sources.append(source)
        new_statuses.append(status)
        records.extend(collected)
        action_status = status.status if status.status in {"ok", "partial"} else "error"
        items.append(
            CompletionItem(
                str(database), analysis.gene, source, "queried", action_status,
                status.record_count, status.message,
            )
        )
        _save_completion_status(
            database, analysis.gene, source, identity, action_status,
            status.record_count, status.message,
        )
    if not queried_sources:
        return analysis, items
    merged_records = _deduplicate_evidence(records)
    merged_statuses = _merge_statuses(analysis.source_statuses, new_statuses)
    enriched = replace(
        analysis,
        generated_at=_utc_now(),
        evidence_record_count=len(merged_records),
        evidence_records=merged_records,
        source_statuses=merged_statuses,
        sources=tuple(
            dict.fromkeys(
                (*analysis.sources, *(status.url for status in new_statuses if status.url))
            )
        ),
    )
    with ExecutionDatabase(database) as execution_database:
        execution_database.save(
            ANALYSIS_KIND, analysis.gene, base_parameters, enriched
        )
    return enriched, items


def _complete_therapy(
    database: Path,
    analysis: AnalysisResult,
    snapshots: tuple[_TherapySnapshot, ...],
    selected: tuple[str, ...],
    *,
    oncokb_token: str | None,
    regulatory_files: tuple[tuple[str, str | Path], ...],
    max_records: int,
    offline: bool,
    dry_run: bool,
    retry_errors: bool,
    strict: bool,
    injected: Sequence[TherapyProvider] | None,
) -> list[CompletionItem]:
    """Role: Backfill missing therapy data for one stored gene."""
    custom = {provider.name.lower(): provider for provider in injected or ()}
    requested = [source for source in selected if source in THERAPY_SOURCES or source in custom]
    requested.extend(source.strip().lower() for source, _path in regulatory_files)
    requested = list(dict.fromkeys(requested))
    if not requested:
        return []
    existing_result = _merge_matching_therapy(analysis, snapshots)
    statuses = _latest_statuses(existing_result.source_statuses)
    records = list(existing_result.associations)
    new_statuses: list[SourceStatus] = []
    queried_sources: list[str] = []
    items: list[CompletionItem] = []
    regulatory_map = {source.strip().lower(): Path(path) for source, path in regulatory_files}
    for source in requested:
        existing = statuses.get(source)
        if existing is not None and (
            _status_is_complete(existing)
            or (existing.status == "error" and not retry_errors)
        ):
            items.append(
                CompletionItem(
                    str(database), analysis.gene, source, "cache_hit",
                    existing.status, existing.record_count, existing.message,
                )
            )
            continue
        if dry_run or offline:
            items.append(
                CompletionItem(
                    str(database), analysis.gene, source,
                    "would_query" if dry_run else "offline_missing", "missing",
                    message="web/local query skipped",
                )
            )
            continue
        provider = custom.get(source) or _therapy_provider(
            source,
            max_records=max_records,
            oncokb_token=oncokb_token,
            regulatory_path=regulatory_map.get(source),
        )
        identity = _provider_identity(provider)
        prior = _load_completion_status(database, analysis.gene, source, identity)
        if prior and (
            prior[0] in COMPLETE_STATUSES
            or (prior[0] == "error" and not retry_errors)
        ):
            items.append(
                CompletionItem(str(database), analysis.gene, source, "cache_hit", *prior)
            )
            continue
        collected, batch_statuses = collect_therapy_associations(
            TherapyContext.from_analysis(analysis), (provider,), strict=strict
        )
        status = batch_statuses[0]
        queried_sources.append(source)
        new_statuses.append(status)
        records.extend(collected)
        action_status = status.status if status.status in {"ok", "partial"} else "error"
        items.append(
            CompletionItem(
                str(database), analysis.gene, source, "queried", action_status,
                status.record_count, status.message,
            )
        )
        _save_completion_status(
            database, analysis.gene, source, identity, action_status,
            status.record_count, status.message,
        )
    if queried_sources:
        merged_records = _deduplicate_therapy(records)
        merged_statuses = _merge_statuses(existing_result.source_statuses, new_statuses)
        result = TherapyResult(
            gene=analysis.gene,
            generated_at=_utc_now(),
            protein_changes=tuple(item.protein_change for item in analysis.shared_mutations),
            variant_therapy_evidence=tuple(
                item for item in merged_records if item.layer == VARIANT_THERAPY_EVIDENCE
            ),
            regulatory_therapy=tuple(
                item for item in merged_records if item.layer == REGULATORY_THERAPY
            ),
            trial_and_target_discovery=tuple(
                item for item in merged_records if item.layer == TRIAL_AND_TARGET_DISCOVERY
            ),
            source_statuses=merged_statuses,
        )
        parameters = next(
            (
                snapshot.parameters
                for snapshot in snapshots
                if _therapy_signature(snapshot.parameters) == _analysis_signature(analysis)
            ),
            {
                "arguments": {
                    "analysis_signature": _analysis_signature(analysis),
                    "completion_sources": sorted(queried_sources),
                }
            },
        )
        with ExecutionDatabase(database) as execution_database:
            execution_database.save(THERAPY_KIND, analysis.gene, parameters, result)
    return items


def _complete_string(
    database: Path,
    gene: str,
    *,
    required_score: int,
    added_nodes: int,
    network_type: str,
    offline: bool,
    dry_run: bool,
    retry_errors: bool,
    http: StringHttp | None,
) -> CompletionItem:
    """Role: Backfill missing string data for one stored gene."""
    source = "string"
    identity = {
        "version": "12.0", "score": required_score,
        "partners": added_nodes, "network_type": network_type,
    }
    cached = _stored_string_counts(
        database, gene, required_score, network_type, added_nodes
    )
    if cached is not None:
        return CompletionItem(
            str(database), gene, source, "cache_hit", "ok",
            cached[1], f"{cached[0]} nodes",
        )
    prior = _load_completion_status(database, gene, source, identity)
    if prior and (
        prior[0] in COMPLETE_STATUSES or (prior[0] == "error" and not retry_errors)
    ):
        return CompletionItem(str(database), gene, source, "cache_hit", *prior)
    if dry_run or offline:
        return CompletionItem(
            str(database), gene, source,
            "would_query" if dry_run else "offline_missing", "missing",
            message="STRING network is missing",
        )
    try:
        run = get_string_network(
            gene, database=database, required_score=required_score,
            network_type=network_type, added_nodes=added_nodes, http=http,
        )
    except (StringNetworkError, OSError, ValueError) as exc:
        _save_completion_status(database, gene, source, identity, "error", 0, str(exc))
        return CompletionItem(str(database), gene, source, "queried", "error", message=str(exc))
    message = f"{len(run.network.nodes)} nodes"
    _save_completion_status(
        database, gene, source, identity, "ok", len(run.network.edges), message
    )
    return CompletionItem(
        str(database), gene, source, "queried", "ok", len(run.network.edges), message
    )


def _complete_uniprot(
    database: Path,
    gene: str,
    *,
    offline: bool,
    dry_run: bool,
    retry_errors: bool,
    client: ProteinClient | None,
) -> CompletionItem:
    """Role: Backfill missing uniprot data for one stored gene."""
    source = "uniprot"
    identity = {"dataset": "reviewed human canonical protein sequence"}
    cached = _stored_sequence(database, gene)
    if cached is not None:
        return CompletionItem(
            str(database), gene, source, "cache_hit", "ok", 1,
            f"{len(cached[0])} residues",
        )
    prior = _load_completion_status(database, gene, source, identity)
    if prior and (
        prior[0] in COMPLETE_STATUSES or (prior[0] == "error" and not retry_errors)
    ):
        return CompletionItem(str(database), gene, source, "cache_hit", *prior)
    if dry_run or offline:
        return CompletionItem(
            str(database), gene, source,
            "would_query" if dry_run else "offline_missing", "missing",
            message="UniProt sequence is missing",
        )
    try:
        annotation = (client or UniProtProteinClient()).fetch(gene)
        if annotation is None or not annotation.sequence:
            _save_completion_status(
                database, gene, source, identity, "unavailable", 0,
                "UniProt did not return a reviewed human entry with a sequence",
            )
            return CompletionItem(
                str(database), gene, source, "queried", "unavailable", 0,
                "UniProt did not return a reviewed human entry with a sequence",
            )
        ProteinSequenceCache(database).store(annotation)
    except (ApiError, OSError, ValueError) as exc:
        _save_completion_status(database, gene, source, identity, "error", 0, str(exc))
        return CompletionItem(str(database), gene, source, "queried", "error", message=str(exc))
    message = f"{len(annotation.sequence)} residues"
    _save_completion_status(database, gene, source, identity, "ok", 1, message)
    return CompletionItem(str(database), gene, source, "queried", "ok", 1, message)


def _evidence_provider(
    source: str, database: Path, gnomad_max_variants: int
) -> EvidenceProvider:
    """Role: Compute evidence provider for the completion workflow."""
    if source == "alphamissense":
        return AlphaMissenseEvidenceProvider(database=database)
    if source == "gnomad":
        return GnomadEvidenceProvider(max_variants=gnomad_max_variants)
    if source == "panelapp":
        return PanelAppEvidenceProvider()
    raise ValueError(f"Unknown evidence provider: {source}")


def _therapy_provider(
    source: str,
    *,
    max_records: int,
    oncokb_token: str | None,
    regulatory_path: Path | None,
) -> TherapyProvider:
    """Role: Compute therapy provider for the completion workflow."""
    if source == "civic":
        return CivicTherapyProvider(max_variants=max_records)
    if source == "oncokb" and oncokb_token:
        return OncoKBTherapyProvider(oncokb_token, max_variants=max_records)
    if source == "clinicaltrials":
        return ClinicalTrialsTherapyProvider(max_studies=max_records)
    if source == "opentargets":
        return OpenTargetsTherapyProvider(max_drugs=max_records)
    if regulatory_path is not None:
        return RegulatoryTableTherapyProvider(
            source,
            regulatory_path,
            jurisdiction={"fda": "US", "ema": "EU"}.get(source),
        )
    raise ValueError(f"Unknown therapy provider: {source}")


def _evidence_context(result: AnalysisResult) -> EvidenceContext:
    """Role: Compute evidence context for the completion workflow."""
    reverse_accessions = {chromosome: accession for accession, chromosome in GRCH38_ACCESSIONS.items()}
    rare: dict[tuple[str, str], RareVariant] = {}
    cancer: dict[tuple[str, str, str], CancerObservation] = {}
    for record in result.evidence_records:
        if not _exact_identity(record):
            continue
        protein_change = record.protein_change or ""
        chromosome = str(record.chromosome)
        accession = reverse_accessions.get(chromosome)
        if not accession or not protein_change:
            continue
        reference = str(record.reference)
        alternate = str(record.alternate)
        position = int(record.position or 0)
        if record.source.casefold() == "clinvar":
            spdi = f"{accession}:{position - 1}:{reference}:{alternate}"
            key = (spdi, protein_change)
            details = dict(record.details)
            disease_name = record.disease_name or ""
            orphanet = (
                record.disease_id.split(":", 1)[1]
                if record.disease_id and record.disease_id.upper().startswith("ORPHA:")
                else ""
            )
            rare.setdefault(
                key,
                RareVariant(
                    gene=result.gene,
                    protein_change=protein_change,
                    raw_protein_change=protein_change,
                    clinvar_id=details.get("clinvar_id", record.source_record_id),
                    clinvar_accession=record.source_record_id.split("|", 1)[0],
                    diseases=(disease_name,) if disease_name else (),
                    orphanet_codes=(orphanet,) if orphanet else (),
                    condition_orphanet_pairs=(
                        ((disease_name, orphanet),) if disease_name and orphanet else ()
                    ),
                    germline_classification=record.classification or "",
                    molecular_consequences=(),
                    canonical_spdi=spdi,
                ),
            )
        elif record.source.casefold() == "gdc":
            details = dict(record.details)
            change = f"chr{chromosome}:g.{position}{reference}>{alternate}"
            case_id = details.get("case_id", record.source_record_id)
            ssm_id = details.get("ssm_id", record.source_record_id.split("|", 1)[0])
            key = (ssm_id, case_id, protein_change)
            cancer.setdefault(
                key,
                CancerObservation(
                    gene=result.gene,
                    protein_change=protein_change,
                    raw_protein_change=protein_change,
                    ssm_id=ssm_id,
                    case_id=case_id,
                    patient_id=details.get("patient_id", case_id),
                    project_id=details.get("project_id", ""),
                    transcript_id=record.transcript,
                    genomic_dna_change=change,
                ),
            )
    for shared in result.shared_mutations:
        for spdi in shared.canonical_spdis:
            rare.setdefault(
                (spdi, shared.protein_change),
                RareVariant(
                    gene=result.gene,
                    protein_change=shared.protein_change,
                    raw_protein_change=shared.protein_change,
                    clinvar_id=shared.clinvar_ids[0] if shared.clinvar_ids else "",
                    clinvar_accession=(
                        shared.clinvar_accessions[0] if shared.clinvar_accessions else ""
                    ),
                    diseases=shared.diseases,
                    orphanet_codes=shared.orphanet_codes,
                    condition_orphanet_pairs=(),
                    germline_classification=(
                        shared.germline_classifications[0]
                        if shared.germline_classifications else ""
                    ),
                    molecular_consequences=(),
                    canonical_spdi=spdi,
                ),
            )
    return EvidenceContext(result.gene, tuple(rare.values()), tuple(cancer.values()))


def _exact_identity(record: EvidenceRecord) -> bool:
    """Role: Compute exact identity for the completion workflow."""
    return bool(
        record.assembly == "GRCh38"
        and record.chromosome
        and record.position
        and record.reference
        and record.alternate
        and set(record.reference.upper()) <= set("ACGT")
        and set(record.alternate.upper()) <= set("ACGT")
    )


def _merge_matching_therapy(
    analysis: AnalysisResult, snapshots: tuple[_TherapySnapshot, ...]
) -> TherapyResult:
    """Role: Merge matching therapy records without losing source provenance."""
    signature = _analysis_signature(analysis)
    chosen_statuses: dict[str, SourceStatus] = {}
    chosen_records: list[TherapyAssociation] = []
    for snapshot in snapshots:
        if _therapy_signature(snapshot.parameters) != signature:
            continue
        sources_here = {status.source.casefold() for status in snapshot.result.source_statuses}
        newly_selected = sources_here - chosen_statuses.keys()
        if not newly_selected:
            continue
        for status in snapshot.result.source_statuses:
            key = status.source.casefold()
            if key in newly_selected:
                chosen_statuses[key] = status
        chosen_records.extend(
            record
            for record in snapshot.result.associations
            if record.source.casefold() in newly_selected
        )
    records = _deduplicate_therapy(chosen_records)
    return TherapyResult(
        gene=analysis.gene,
        generated_at=_utc_now(),
        protein_changes=tuple(item.protein_change for item in analysis.shared_mutations),
        variant_therapy_evidence=tuple(
            item for item in records if item.layer == VARIANT_THERAPY_EVIDENCE
        ),
        regulatory_therapy=tuple(
            item for item in records if item.layer == REGULATORY_THERAPY
        ),
        trial_and_target_discovery=tuple(
            item for item in records if item.layer == TRIAL_AND_TARGET_DISCOVERY
        ),
        source_statuses=tuple(chosen_statuses.values()),
    )


def _therapy_signature(parameters: Mapping[str, Any]) -> Mapping[str, Any] | None:
    """Role: Compute therapy signature for the completion workflow."""
    arguments = parameters.get("arguments")
    if isinstance(arguments, Mapping):
        signature = arguments.get("analysis_signature")
        if isinstance(signature, Mapping):
            return signature
    signature = parameters.get("analysis_signature")
    return signature if isinstance(signature, Mapping) else None


def _analysis_signature(result: AnalysisResult) -> dict[str, Any]:
    """Role: Compute analysis signature for the completion workflow."""
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


def _latest_statuses(statuses: Iterable[SourceStatus]) -> dict[str, SourceStatus]:
    """Role: Compute latest statuses for the completion workflow."""
    output: dict[str, SourceStatus] = {}
    for status in reversed(tuple(statuses)):
        output.setdefault(status.source.casefold(), status)
    return output


def _status_is_complete(status: SourceStatus | None) -> bool:
    """Role: Compute status is complete for the completion workflow."""
    return status is not None and status.status in {"ok", "partial"}


def _merge_statuses(
    existing: Iterable[SourceStatus], new: Iterable[SourceStatus]
) -> tuple[SourceStatus, ...]:
    """Role: Merge statuses records without losing source provenance."""
    merged = _latest_statuses(existing)
    for status in new:
        merged[status.source.casefold()] = status
    return tuple(sorted(merged.values(), key=lambda item: item.source.casefold()))


def _deduplicate_evidence(records: Iterable[EvidenceRecord]) -> tuple[EvidenceRecord, ...]:
    """Role: Deduplicate evidence while preserving provenance."""
    unique = {
        (
            item.source, item.source_record_id, item.evidence_layer,
            item.variant_key, item.disease_id,
        ): item
        for item in records
    }
    return tuple(
        sorted(
            unique.values(),
            key=lambda item: (
                item.evidence_layer, item.source, item.variant_key or "",
                item.source_record_id,
            ),
        )
    )


def _deduplicate_therapy(
    records: Iterable[TherapyAssociation],
) -> tuple[TherapyAssociation, ...]:
    """Role: Deduplicate therapy while preserving provenance."""
    unique = {
        (
            item.layer, item.source, item.source_record_id, item.therapy_name,
            item.protein_change, item.molecular_profile, item.disease_ids,
            item.disease_names, item.effect, item.therapy_status, item.jurisdiction,
        ): item
        for item in records
    }
    return tuple(
        sorted(
            unique.values(),
            key=lambda item: (
                item.layer, item.source, item.protein_change or "",
                item.therapy_name.casefold(), item.source_record_id,
            ),
        )
    )


def _initialize_completion_status(database: Path) -> None:
    """Role: Compute initialize completion status for the completion workflow."""
    with sqlite3.connect(database) as connection:
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS completion_query_status (
                gene TEXT NOT NULL,
                source TEXT NOT NULL,
                request_hash TEXT NOT NULL,
                status TEXT NOT NULL,
                record_count INTEGER NOT NULL,
                message TEXT,
                queried_at TEXT NOT NULL,
                schema_version INTEGER NOT NULL,
                PRIMARY KEY (gene, source, request_hash)
            )
            """
        )
    _restrict_permissions(database)


def _stored_string_counts(
    database: Path,
    gene: str,
    required_score: int,
    network_type: str,
    added_nodes: int,
) -> tuple[int, int] | None:
    """Role: Compute stored string counts for the completion workflow."""
    try:
        with sqlite3.connect(f"{database.as_uri()}?mode=ro", uri=True) as connection:
            exists = connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='string_network_cache'"
            ).fetchone()
            if not exists:
                return None
            row = connection.execute(
                """
                SELECT result_json FROM string_network_cache
                WHERE gene = ? AND species = 9606 AND string_version = '12.0'
                  AND required_score = ? AND network_type = ? AND added_nodes = ?
                """,
                (gene.upper(), required_score, network_type, added_nodes),
            ).fetchone()
    except sqlite3.Error:
        return None
    if row is None:
        return None
    try:
        payload = json.loads(str(row[0]))
    except json.JSONDecodeError:
        return None
    if not isinstance(payload, Mapping):
        return None
    nodes = payload.get("nodes")
    edges = payload.get("edges")
    return (
        len(nodes) if isinstance(nodes, list) else 0,
        len(edges) if isinstance(edges, list) else 0,
    )


def _stored_sequence(database: Path, gene: str) -> tuple[str, str | None] | None:
    """Role: Compute stored sequence for the completion workflow."""
    try:
        with sqlite3.connect(f"{database.as_uri()}?mode=ro", uri=True) as connection:
            exists = connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='protein_sequence_cache'"
            ).fetchone()
            if not exists:
                return None
            row = connection.execute(
                "SELECT sequence, accession FROM protein_sequence_cache WHERE gene = ?",
                (gene.upper(),),
            ).fetchone()
    except sqlite3.Error:
        return None
    if row is None:
        return None
    sequence = str(row[0]).strip().upper()
    if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        return None
    return sequence, str(row[1]) if row[1] else None


def _load_completion_status(
    database: Path,
    gene: str,
    source: str,
    identity: Mapping[str, Any],
) -> tuple[str, int, str | None] | None:
    """Role: Compute load completion status for the completion workflow."""
    request_hash = _hash_identity(identity)
    try:
        with sqlite3.connect(f"{database.as_uri()}?mode=ro", uri=True) as connection:
            exists = connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='completion_query_status'"
            ).fetchone()
            if not exists:
                return None
            row = connection.execute(
                """
                SELECT status, record_count, message
                FROM completion_query_status
                WHERE gene = ? AND source = ? AND request_hash = ?
                  AND schema_version = ?
                """,
                (gene.upper(), source.lower(), request_hash, COMPLETION_SCHEMA_VERSION),
            ).fetchone()
    except sqlite3.Error:
        return None
    if row is None:
        return None
    return str(row[0]), int(row[1]), str(row[2]) if row[2] is not None else None


def _save_completion_status(
    database: Path,
    gene: str,
    source: str,
    identity: Mapping[str, Any],
    status: str,
    record_count: int,
    message: str | None,
) -> None:
    """Role: Compute save completion status for the completion workflow."""
    _initialize_completion_status(database)
    with sqlite3.connect(database) as connection:
        connection.execute(
            """
            INSERT INTO completion_query_status (
                gene, source, request_hash, status, record_count, message,
                queried_at, schema_version
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(gene, source, request_hash) DO UPDATE SET
                status = excluded.status,
                record_count = excluded.record_count,
                message = excluded.message,
                queried_at = excluded.queried_at,
                schema_version = excluded.schema_version
            """,
            (
                gene.upper(), source.lower(), _hash_identity(identity), status,
                record_count, message, _utc_now(), COMPLETION_SCHEMA_VERSION,
            ),
        )
    _restrict_permissions(database)


def _provider_identity(provider: Any) -> dict[str, Any]:
    """Role: Resolve or describe the configured identity provider."""
    method = getattr(provider, "cache_identity", None)
    if not callable(method):
        raise TypeError(
            f"{type(provider).__name__} must implement cache_identity() for backfill"
        )
    identity = method()
    if not isinstance(identity, Mapping):
        raise TypeError("cache_identity() must return a mapping")
    return {
        "type": f"{type(provider).__module__}.{type(provider).__qualname__}",
        "configuration": _redact_identity(identity),
    }


def _redact_identity(value: Any) -> Any:
    """Role: Compute redact identity for the completion workflow."""
    if isinstance(value, Mapping):
        output: dict[str, Any] = {}
        for key, item in value.items():
            folded = str(key).casefold()
            output[str(key)] = (
                "REDACTED"
                if any(marker in folded for marker in ("token", "secret", "password", "credential"))
                else _redact_identity(item)
            )
        return output
    if isinstance(value, Path):
        resolved = value.expanduser().resolve()
        descriptor: dict[str, Any] = {"path": str(resolved)}
        if resolved.is_file():
            descriptor["sha256"] = _file_sha256(resolved)
        return descriptor
    if isinstance(value, (tuple, list)):
        return [_redact_identity(item) for item in value]
    return value


def _hash_identity(identity: Mapping[str, Any]) -> str:
    """Role: Compute a stable hash for identity."""
    payload = json.dumps(
        _redact_identity(identity), ensure_ascii=False, sort_keys=True,
        separators=(",", ":"), allow_nan=False,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _file_sha256(path: Path) -> str:
    """Role: Compute or inspect sha256 file metadata."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_csv(value: Any) -> Any:
    """Role: Sanitize csv for safe serialization."""
    if value is None or not isinstance(value, str):
        return value
    if value.startswith(("\t", "\r", "\n")) or value.lstrip(" \t\r\n").startswith(
        ("=", "+", "-", "@")
    ):
        return "'" + value
    return value


def _restrict_permissions(path: Path) -> None:
    """Role: Compute restrict permissions for the completion workflow."""
    if os.name != "posix":
        return
    for candidate in (path, Path(f"{path}-wal"), Path(f"{path}-shm")):
        try:
            os.chmod(candidate, 0o600)
        except FileNotFoundError:
            continue


def _utc_now() -> str:
    """Role: Return the current UTC timestamp in canonical form."""
    return datetime.now(timezone.utc).isoformat()
