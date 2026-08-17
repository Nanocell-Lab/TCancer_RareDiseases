"""Incremental cross-gene analysis focused on actionable protein families."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import sqlite3
from collections import Counter
from contextlib import closing
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping


GLOBAL_INDEX_VERSION = 1
ACTIVE_REGULATORY_STATUSES = {
    "regulatory_approved",
    "label_biomarker",
    "conditional_approval",
}


@dataclass(frozen=True)
class ProteinFamilyAssignment:
    gene: str
    family_group: str
    protein_family: str


@dataclass(frozen=True)
class GlobalGeneSummary:
    gene: str
    family_group: str
    protein_family: str
    analysis_snapshot_key: str
    therapy_snapshot_key: str | None
    rare_variant_count: int
    cancer_observation_count: int
    shared_variant_count: int
    confirmed_shared_count: int
    mixed_shared_count: int
    candidate_shared_count: int
    max_variant_patient_count: int
    rare_disease_count: int
    cancer_project_count: int
    variant_specific_therapy_count: int
    active_regulatory_therapy_count: int
    clinical_trial_count: int
    mechanism_target_count: int
    actionability_tier: str
    source_coverage_count: int


@dataclass(frozen=True)
class GlobalFamilySummary:
    family_group: str
    protein_family: str
    gene_count: int
    genes_with_shared_variants: int
    shared_variant_count: int
    confirmed_shared_count: int
    genes_with_variant_therapy: int
    genes_with_active_regulatory: int
    genes_with_clinical_trials: int
    genes_with_mechanism_target: int
    rare_disease_count: int
    cancer_project_count: int


@dataclass(frozen=True)
class GlobalFamilyGroupSummary:
    family_group: str
    gene_count: int
    genes_with_shared_variants: int
    shared_variant_count: int
    confirmed_shared_count: int
    genes_with_variant_therapy: int
    genes_with_active_regulatory: int
    genes_with_clinical_trials: int
    genes_with_mechanism_target: int


@dataclass(frozen=True)
class FigureProposal:
    figure_id: str
    title: str
    purpose: str
    implemented: bool
    interaction: str


@dataclass(frozen=True)
class GlobalFamilyReport:
    generated_at: str
    index_database: str
    source_databases: tuple[str, ...]
    indexed_new_snapshots: int
    active_analysis_count: int
    active_therapy_count: int
    unmatched_therapy_gene_count: int
    index_error_count: int
    source_error_count: int
    gene_summaries: tuple[GlobalGeneSummary, ...]
    family_group_summaries: tuple[GlobalFamilyGroupSummary, ...]
    family_summaries: tuple[GlobalFamilySummary, ...]
    disease_family_counts: tuple[tuple[str, str, int], ...]
    project_family_counts: tuple[tuple[str, str, int], ...]
    figure_proposals: tuple[FigureProposal, ...]
    caveat: str = (
        "Actionability remains layered: variant-specific evidence, regulatory "
        "status, clinical trials, and mechanism/target evidence are not equivalent."
    )

    def to_dict(self) -> dict[str, Any]:
        """Role: Convert the domain object into a JSON-safe dictionary."""
        return asdict(self)


class GlobalAnalysisIndex:
    """Normalize each new execution snapshot once and reuse compact projections."""

    def __init__(self, path: str | Path) -> None:
        """Role: Compute init for the global_analysis workflow."""
        self.path = Path(path).expanduser().resolve()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if os.name == "posix" and not self.path.exists():
            descriptor = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_RDWR, 0o600)
            os.close(descriptor)
        self._initialize()

    def sync(self, databases: Iterable[str | Path]) -> int:
        """Role: Provide the sync operation for the global_analysis workflow."""
        indexed = 0
        for value in _unique_paths(databases):
            if value == self.path or not value.is_file():
                continue
            indexed += self._sync_database(value)
        self._restrict_permissions()
        return indexed

    def build_report(
        self,
        databases: Iterable[str | Path],
        family_catalog: Mapping[str, ProteinFamilyAssignment],
        *,
        indexed_new_snapshots: int = 0,
    ) -> GlobalFamilyReport:
        """Role: Build report data for downstream use."""
        source_paths = tuple(str(path) for path in _unique_paths(databases) if path != self.path)
        with closing(sqlite3.connect(self.path)) as connection:
            connection.row_factory = sqlite3.Row
            analyses = connection.execute(
                """
                SELECT * FROM global_snapshot
                WHERE kind = 'mutation_analysis'
                ORDER BY gene, created_at DESC, source_execution_id DESC
                """
            ).fetchall()
            therapies = connection.execute(
                """
                SELECT * FROM global_snapshot
                WHERE kind = 'therapy_analysis'
                ORDER BY gene, created_at DESC, source_execution_id DESC
                """
            ).fetchall()
            index_error_count = int(
                connection.execute("SELECT COUNT(*) FROM global_index_error").fetchone()[0]
            )
            source_error_count = int(
                connection.execute(
                    "SELECT COUNT(*) FROM global_source_cursor WHERE last_error IS NOT NULL"
                ).fetchone()[0]
            )
        active_analyses = _latest_by_gene(analyses)
        therapies_by_gene: dict[str, list[sqlite3.Row]] = {}
        for row in therapies:
            therapies_by_gene.setdefault(str(row["gene"]), []).append(row)

        gene_summaries: list[GlobalGeneSummary] = []
        disease_family: Counter[tuple[str, str]] = Counter()
        project_family: Counter[tuple[str, str]] = Counter()
        unmatched_therapy = 0
        active_therapy_count = 0
        for gene, analysis_row in sorted(active_analyses.items()):
            analysis = _json_object(analysis_row["summary_json"])
            matching_therapy = next(
                (
                    row
                    for row in therapies_by_gene.get(gene, [])
                    if row["therapy_signature_hash"]
                    and row["therapy_signature_hash"]
                    == analysis_row["analysis_signature_hash"]
                ),
                None,
            )
            therapy: Mapping[str, Any] = {}
            therapy_key: str | None = None
            if matching_therapy is not None:
                therapy = _json_object(matching_therapy["summary_json"])
                therapy_key = str(matching_therapy["snapshot_key"])
                active_therapy_count += 1
            elif therapies_by_gene.get(gene):
                unmatched_therapy += 1
            assignment = family_catalog.get(
                gene,
                ProteinFamilyAssignment(
                    gene, "Other analyzed proteins", "Unclassified"
                ),
            )
            for disease, count in _pairs(analysis.get("disease_counts")):
                disease_family[(str(disease), assignment.protein_family)] += int(count)
            for project, count in _pairs(analysis.get("project_counts")):
                project_family[(str(project), assignment.protein_family)] += int(count)
            variant_therapy = _integer(therapy.get("variant_specific_count"))
            active_regulatory = _integer(therapy.get("active_regulatory_count"))
            clinical_trials = _integer(therapy.get("clinical_trial_count"))
            mechanism_targets = _integer(therapy.get("mechanism_target_count"))
            gene_summaries.append(
                GlobalGeneSummary(
                    gene=gene,
                    family_group=assignment.family_group,
                    protein_family=assignment.protein_family,
                    analysis_snapshot_key=str(analysis_row["snapshot_key"]),
                    therapy_snapshot_key=therapy_key,
                    rare_variant_count=_integer(analysis.get("rare_variant_count")),
                    cancer_observation_count=_integer(
                        analysis.get("cancer_observation_count")
                    ),
                    shared_variant_count=_integer(analysis.get("shared_variant_count")),
                    confirmed_shared_count=_integer(analysis.get("confirmed_count")),
                    mixed_shared_count=_integer(analysis.get("mixed_count")),
                    candidate_shared_count=_integer(analysis.get("candidate_count")),
                    max_variant_patient_count=_integer(
                        analysis.get("max_variant_patient_count")
                    ),
                    rare_disease_count=len(_pairs(analysis.get("disease_counts"))),
                    cancer_project_count=len(_pairs(analysis.get("project_counts"))),
                    variant_specific_therapy_count=variant_therapy,
                    active_regulatory_therapy_count=active_regulatory,
                    clinical_trial_count=clinical_trials,
                    mechanism_target_count=mechanism_targets,
                    actionability_tier=_actionability_tier(
                        variant_therapy,
                        active_regulatory,
                        clinical_trials,
                        mechanism_targets,
                    ),
                    source_coverage_count=_integer(
                        analysis.get("source_coverage_count")
                    ),
                )
            )
        families = _family_summaries(gene_summaries, disease_family, project_family)
        return GlobalFamilyReport(
            generated_at=datetime.now(timezone.utc).isoformat(),
            index_database=str(self.path),
            source_databases=source_paths,
            indexed_new_snapshots=indexed_new_snapshots,
            active_analysis_count=len(active_analyses),
            active_therapy_count=active_therapy_count,
            unmatched_therapy_gene_count=unmatched_therapy,
            index_error_count=index_error_count,
            source_error_count=source_error_count,
            gene_summaries=tuple(gene_summaries),
            family_group_summaries=tuple(_family_group_summaries(gene_summaries)),
            family_summaries=tuple(families),
            disease_family_counts=tuple(
                (disease, family, count)
                for (disease, family), count in sorted(
                    disease_family.items(), key=lambda item: (-item[1], item[0])
                )
            ),
            project_family_counts=tuple(
                (project, family, count)
                for (project, family), count in sorted(
                    project_family.items(), key=lambda item: (-item[1], item[0])
                )
            ),
            figure_proposals=_figure_proposals(),
        )

    def _sync_database(self, source: Path) -> int:
        """Role: Compute sync database for the global_analysis workflow."""
        source_text = str(source)
        with closing(sqlite3.connect(self.path)) as index:
            cursor_row = index.execute(
                "SELECT last_execution_id FROM global_source_cursor WHERE source_path = ?",
                (source_text,),
            ).fetchone()
            last_id = int(cursor_row[0]) if cursor_row else 0
        try:
            with closing(_connect_source(source)) as connection:
                if not _has_execution_cache(connection):
                    return 0
                maximum_row = connection.execute(
                    "SELECT COALESCE(MAX(id), 0) FROM execution_cache"
                ).fetchone()
                maximum_id = int(maximum_row[0])
                if maximum_id < last_id:
                    last_id = 0
                if maximum_id == last_id:
                    self._touch_source(source_text, last_id, None)
                    return 0
                rows = connection.execute(
                    """
                    SELECT id, kind, gene, parameter_hash, parameters_json,
                           result_json, created_at
                    FROM execution_cache
                    WHERE id > ?
                    ORDER BY id
                    """,
                    (last_id,),
                ).fetchall()
        except sqlite3.Error as exc:
            self._touch_source(source_text, last_id, str(exc))
            return 0

        indexed = 0
        with closing(sqlite3.connect(self.path)) as index:
            with index:
                for row in rows:
                    try:
                        result = _json_object(row["result_json"])
                        parameters = _json_object(row["parameters_json"])
                        kind = str(row["kind"])
                        summary, analysis_signature, therapy_signature = _project_snapshot(
                            kind, result, parameters
                        )
                    except (TypeError, ValueError, json.JSONDecodeError) as exc:
                        index.execute(
                            """
                            INSERT INTO global_index_error (
                                source_path, source_execution_id, message, observed_at
                            ) VALUES (?, ?, ?, ?)
                            """,
                            (source_text, int(row["id"]), str(exc), _utc_now()),
                        )
                        continue
                    result_hash = hashlib.sha256(
                        str(row["result_json"]).encode("utf-8")
                    ).hexdigest()
                    snapshot_key = _hash_object(
                        {
                            "kind": kind,
                            "gene": str(row["gene"]).upper(),
                            "parameter_hash": str(row["parameter_hash"]),
                            "result_hash": result_hash,
                            "created_at": str(row["created_at"]),
                        }
                    )
                    before = index.total_changes
                    index.execute(
                        """
                        INSERT OR IGNORE INTO global_snapshot (
                            snapshot_key, source_path, source_execution_id, kind,
                            gene, parameter_hash, result_hash, created_at,
                            analysis_signature_hash, therapy_signature_hash,
                            summary_json, projection_version, indexed_at
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            snapshot_key,
                            source_text,
                            int(row["id"]),
                            kind,
                            str(row["gene"]).upper(),
                            str(row["parameter_hash"]),
                            result_hash,
                            str(row["created_at"]),
                            analysis_signature,
                            therapy_signature,
                            _canonical_json(summary),
                            GLOBAL_INDEX_VERSION,
                            _utc_now(),
                        ),
                    )
                    indexed += int(index.total_changes > before)
                index.execute(
                    """
                    INSERT INTO global_source_cursor (
                        source_path, last_execution_id, last_scanned_at, last_error
                    ) VALUES (?, ?, ?, NULL)
                    ON CONFLICT(source_path) DO UPDATE SET
                        last_execution_id = excluded.last_execution_id,
                        last_scanned_at = excluded.last_scanned_at,
                        last_error = NULL
                    """,
                    (source_text, maximum_id, _utc_now()),
                )
        return indexed

    def _touch_source(self, source: str, last_id: int, error: str | None) -> None:
        """Role: Compute touch source for the global_analysis workflow."""
        with closing(sqlite3.connect(self.path)) as connection:
            with connection:
                connection.execute(
                    """
                    INSERT INTO global_source_cursor (
                        source_path, last_execution_id, last_scanned_at, last_error
                    ) VALUES (?, ?, ?, ?)
                    ON CONFLICT(source_path) DO UPDATE SET
                        last_execution_id = excluded.last_execution_id,
                        last_scanned_at = excluded.last_scanned_at,
                        last_error = excluded.last_error
                    """,
                    (source, last_id, _utc_now(), error),
                )

    def _initialize(self) -> None:
        """Role: Compute initialize for the global_analysis workflow."""
        with closing(sqlite3.connect(self.path)) as connection:
            with connection:
                connection.executescript(
                    """
                    PRAGMA journal_mode = WAL;
                    PRAGMA busy_timeout = 5000;
                    CREATE TABLE IF NOT EXISTS global_meta (
                        key TEXT PRIMARY KEY,
                        value TEXT NOT NULL
                    );
                    CREATE TABLE IF NOT EXISTS global_source_cursor (
                        source_path TEXT PRIMARY KEY,
                        last_execution_id INTEGER NOT NULL DEFAULT 0,
                        last_scanned_at TEXT NOT NULL,
                        last_error TEXT
                    );
                    CREATE TABLE IF NOT EXISTS global_snapshot (
                        snapshot_key TEXT PRIMARY KEY,
                        source_path TEXT NOT NULL,
                        source_execution_id INTEGER NOT NULL,
                        kind TEXT NOT NULL,
                        gene TEXT NOT NULL,
                        parameter_hash TEXT NOT NULL,
                        result_hash TEXT NOT NULL,
                        created_at TEXT NOT NULL,
                        analysis_signature_hash TEXT,
                        therapy_signature_hash TEXT,
                        summary_json TEXT NOT NULL,
                        projection_version INTEGER NOT NULL,
                        indexed_at TEXT NOT NULL
                    );
                    CREATE INDEX IF NOT EXISTS idx_global_snapshot_active
                    ON global_snapshot(kind, gene, created_at DESC);
                    CREATE TABLE IF NOT EXISTS global_index_error (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        source_path TEXT NOT NULL,
                        source_execution_id INTEGER NOT NULL,
                        message TEXT NOT NULL,
                        observed_at TEXT NOT NULL
                    );
                    """
                )
                version_row = connection.execute(
                    "SELECT value FROM global_meta WHERE key = 'projection_version'"
                ).fetchone()
                if version_row is not None and int(version_row[0]) != GLOBAL_INDEX_VERSION:
                    connection.execute("DELETE FROM global_snapshot")
                    connection.execute("DELETE FROM global_source_cursor")
                    connection.execute("DELETE FROM global_index_error")
                connection.execute(
                    """
                    INSERT INTO global_meta(key, value) VALUES ('projection_version', ?)
                    ON CONFLICT(key) DO UPDATE SET value = excluded.value
                    """,
                    (str(GLOBAL_INDEX_VERSION),),
                )
        self._restrict_permissions()

    def _restrict_permissions(self) -> None:
        """Role: Compute restrict permissions for the global_analysis workflow."""
        if os.name != "posix":
            return
        for candidate in (self.path, Path(f"{self.path}-wal"), Path(f"{self.path}-shm")):
            try:
                os.chmod(candidate, 0o600)
            except FileNotFoundError:
                continue


def discover_global_databases(
    results_root: str | Path, *, exclude: Iterable[str | Path] = ()
) -> tuple[Path, ...]:
    """Role: Discover global databases resources from configured locations."""
    root = Path(results_root).expanduser().resolve()
    excluded = {Path(value).expanduser().resolve() for value in exclude}
    return tuple(
        sorted(
            (
                path.resolve()
                for path in root.rglob("*.sqlite3")
                if path.is_file() and path.resolve() not in excluded
            ),
            key=str,
        )
    )


def build_global_family_analysis(
    databases: Iterable[str | Path],
    *,
    index_database: str | Path,
    gene_list_dir: str | Path,
) -> GlobalFamilyReport:
    """Role: Build global family analysis data for downstream use."""
    paths = _unique_paths(databases)
    index = GlobalAnalysisIndex(index_database)
    indexed = index.sync(paths)
    return index.build_report(
        paths,
        load_family_catalog(gene_list_dir),
        indexed_new_snapshots=indexed,
    )


def load_family_catalog(path: str | Path) -> dict[str, ProteinFamilyAssignment]:
    """Role: Load and validate family catalog data."""
    root = Path(path)
    groups = (
        ("Receptor tyrosine kinases", "receptor_tyrosine_kinases.txt"),
        ("TGF-beta receptor superfamily", "tgfb_receptors.txt"),
        ("Kinase-linked cytokine receptors", "kinase_linked_cytokine_receptors.txt"),
    )
    catalog: dict[str, ProteinFamilyAssignment] = {}
    for group, filename in groups:
        source = root / filename
        if not source.is_file():
            raise FileNotFoundError(f"Family list does not exist: {source}")
        for gene in source.read_text(encoding="utf-8").splitlines():
            symbol = gene.strip().upper()
            if not symbol:
                continue
            if symbol in catalog:
                raise ValueError(f"Gene {symbol} appears in more than one family")
            catalog[symbol] = ProteinFamilyAssignment(
                symbol, group, _protein_subfamily(symbol, group)
            )
    return catalog


def write_global_report(report: GlobalFamilyReport, output_dir: str | Path) -> tuple[Path, ...]:
    """Role: Serialize and write global report output."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "global_family_analysis.json"
    json_path.write_text(
        json.dumps(report.to_dict(), ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    gene_path = output / "global_gene_summary.csv"
    family_path = output / "global_family_summary.csv"
    group_path = output / "global_family_group_summary.csv"
    _write_dataclass_csv(gene_path, report.gene_summaries)
    _write_dataclass_csv(family_path, report.family_summaries)
    _write_dataclass_csv(group_path, report.family_group_summaries)
    return json_path, gene_path, family_path, group_path


def _project_snapshot(
    kind: str, result: Mapping[str, Any], parameters: Mapping[str, Any]
) -> tuple[dict[str, Any], str | None, str | None]:
    """Role: Compute project snapshot for the global_analysis workflow."""
    if kind == "mutation_analysis":
        shared = _records(result.get("shared_mutations"))
        match_counts = Counter(str(item.get("match_level") or "unknown") for item in shared)
        diseases: Counter[str] = Counter()
        projects: Counter[str] = Counter()
        max_patients = 0
        for item in shared:
            max_patients = max(max_patients, _integer(item.get("cancer_patient_count")))
            for disease in _strings(item.get("diseases")):
                diseases[disease] += 1
            for project in _strings(item.get("tcga_projects")):
                projects[project] += 1
        signature = _analysis_signature_from_result(result)
        statuses = {
            str(item.get("source"))
            for item in _records(result.get("source_statuses"))
            if str(item.get("status")) in {"ok", "partial"}
        }
        return (
            {
                "rare_variant_count": _integer(result.get("rare_variant_count")),
                "cancer_observation_count": _integer(
                    result.get("cancer_observation_count")
                ),
                "shared_variant_count": len(shared),
                "confirmed_count": match_counts["confirmed_genomic_and_protein"],
                "mixed_count": match_counts["mixed_genomic_and_candidate_protein"],
                "candidate_count": match_counts["candidate_protein_match"],
                "max_variant_patient_count": max_patients,
                "disease_counts": sorted(diseases.items()),
                "project_counts": sorted(projects.items()),
                "source_coverage_count": len(statuses),
            },
            _hash_object(signature),
            None,
        )
    if kind == "therapy_analysis":
        variant = _records(result.get("variant_therapy_evidence"))
        regulatory = _records(result.get("regulatory_therapy"))
        discovery = _records(result.get("trial_and_target_discovery"))
        active_regulatory = sum(
            str(item.get("therapy_status")) in ACTIVE_REGULATORY_STATUSES
            for item in regulatory
        )
        clinical_trials = sum(
            str(item.get("therapy_status")) == "clinical_trial" for item in discovery
        )
        mechanism_targets = sum(
            str(item.get("therapy_status")) == "mechanism_target" for item in discovery
        )
        signature = _therapy_signature_from_parameters(parameters)
        return (
            {
                "variant_specific_count": len(variant),
                "active_regulatory_count": active_regulatory,
                "clinical_trial_count": clinical_trials,
                "mechanism_target_count": mechanism_targets,
            },
            None,
            _hash_object(signature) if signature is not None else None,
        )
    raise ValueError(f"Unknown global snapshot type: {kind}")


def _analysis_signature_from_result(result: Mapping[str, Any]) -> dict[str, Any]:
    """Role: Compute analysis signature from result for the global_analysis workflow."""
    return {
        "gene": str(result.get("gene") or "").upper(),
        "shared_mutations": [
            {
                "protein_change": item.get("protein_change"),
                "diseases": list(_strings(item.get("diseases"))),
                "orphanet_codes": list(_strings(item.get("orphanet_codes"))),
                "match_level": item.get("match_level"),
            }
            for item in _records(result.get("shared_mutations"))
        ],
    }


def _therapy_signature_from_parameters(
    parameters: Mapping[str, Any],
) -> Mapping[str, Any] | None:
    """Role: Compute therapy signature from parameters for the global_analysis workflow."""
    arguments = parameters.get("arguments")
    if isinstance(arguments, dict) and isinstance(arguments.get("analysis_signature"), dict):
        return arguments["analysis_signature"]
    if isinstance(parameters.get("analysis_signature"), dict):
        return parameters["analysis_signature"]
    return None


def _family_summaries(
    genes: Iterable[GlobalGeneSummary],
    disease_counts: Counter[tuple[str, str]],
    project_counts: Counter[tuple[str, str]],
) -> list[GlobalFamilySummary]:
    """Role: Compute family summaries for the global_analysis workflow."""
    grouped: dict[tuple[str, str], list[GlobalGeneSummary]] = {}
    for gene in genes:
        grouped.setdefault((gene.family_group, gene.protein_family), []).append(gene)
    output = []
    for (group, family), members in grouped.items():
        output.append(
            GlobalFamilySummary(
                family_group=group,
                protein_family=family,
                gene_count=len(members),
                genes_with_shared_variants=sum(item.shared_variant_count > 0 for item in members),
                shared_variant_count=sum(item.shared_variant_count for item in members),
                confirmed_shared_count=sum(item.confirmed_shared_count for item in members),
                genes_with_variant_therapy=sum(item.variant_specific_therapy_count > 0 for item in members),
                genes_with_active_regulatory=sum(item.active_regulatory_therapy_count > 0 for item in members),
                genes_with_clinical_trials=sum(item.clinical_trial_count > 0 for item in members),
                genes_with_mechanism_target=sum(item.mechanism_target_count > 0 for item in members),
                rare_disease_count=sum(1 for disease, item_family in disease_counts if item_family == family),
                cancer_project_count=sum(1 for project, item_family in project_counts if item_family == family),
            )
        )
    return sorted(
        output,
        key=lambda item: (
            -item.genes_with_variant_therapy,
            -item.genes_with_active_regulatory,
            -item.confirmed_shared_count,
            item.protein_family,
        ),
    )


def _family_group_summaries(
    genes: Iterable[GlobalGeneSummary],
) -> list[GlobalFamilyGroupSummary]:
    """Role: Compute family group summaries for the global_analysis workflow."""
    grouped: dict[str, list[GlobalGeneSummary]] = {}
    for gene in genes:
        grouped.setdefault(gene.family_group, []).append(gene)
    return sorted(
        (
            GlobalFamilyGroupSummary(
                family_group=group,
                gene_count=len(members),
                genes_with_shared_variants=sum(item.shared_variant_count > 0 for item in members),
                shared_variant_count=sum(item.shared_variant_count for item in members),
                confirmed_shared_count=sum(item.confirmed_shared_count for item in members),
                genes_with_variant_therapy=sum(item.variant_specific_therapy_count > 0 for item in members),
                genes_with_active_regulatory=sum(item.active_regulatory_therapy_count > 0 for item in members),
                genes_with_clinical_trials=sum(item.clinical_trial_count > 0 for item in members),
                genes_with_mechanism_target=sum(item.mechanism_target_count > 0 for item in members),
            )
            for group, members in grouped.items()
        ),
        key=lambda item: (-item.genes_with_variant_therapy, -item.confirmed_shared_count, item.family_group),
    )


def _actionability_tier(variant: int, regulatory: int, trials: int, mechanism: int) -> str:
    """Role: Compute actionability tier for the global_analysis workflow."""
    if regulatory:
        return "Active regulatory biomarker"
    if variant:
        return "Variant-specific evidence"
    if trials:
        return "Clinical investigation"
    if mechanism:
        return "Mechanism/target only"
    return "No therapy evidence"


def _figure_proposals() -> tuple[FigureProposal, ...]:
    """Role: Compute figure proposals for the global_analysis workflow."""
    return (
        FigureProposal(
            "family_actionability_heatmap",
            "Protein-family actionability matrix",
            "Compare confirmed rare-oncology overlap and distinct therapy layers by family.",
            True,
            "Filter family group and open the contributing genes.",
        ),
        FigureProposal(
            "gene_family_landscape",
            "Gene-level rare-oncology landscape",
            "Relate confirmed shared variants to maximum recurrence, grouped by protein family.",
            True,
            "Hover gene, family, actionability tier and raw counts.",
        ),
        FigureProposal(
            "family_actionability_stack",
            "Highest actionability tier by family",
            "Show how many genes have variant, regulatory, trial, target-only or no therapy evidence.",
            True,
            "Toggle tiers and select a family for gene-level detail.",
        ),
        FigureProposal(
            "rare_disease_family_heatmap",
            "Rare-disease convergence across protein families",
            "Identify rare diseases whose cancer-shared variants span several receptor families.",
            True,
            "Search disease and limit to confirmed shared variants.",
        ),
        FigureProposal(
            "variant_evidence_upset",
            "Variant evidence intersections",
            "Intersect genomic confirmation, pathogenic classification, population rarity and variant therapy.",
            False,
            "Click an intersection to inspect exact variants and provenance.",
        ),
    )


def _protein_subfamily(gene: str, group: str) -> str:
    """Role: Compute protein subfamily for the global_analysis workflow."""
    if group == "TGF-beta receptor superfamily":
        if gene == "TGFBR3":
            return "TGF-beta accessory receptor"
        if gene in {"ACVR2A", "ACVR2B", "AMHR2", "BMPR2", "TGFBR2"}:
            return "TGF-beta type II receptors"
        return "TGF-beta type I receptors"
    if group == "Kinase-linked cytokine receptors":
        cytokine_groups = (
            ("Common gamma-chain cytokine receptors", {"IL2RA", "IL2RB", "IL2RG", "IL4R", "IL7R", "IL9R", "IL15RA", "IL21R"}),
            ("IL-3 receptor family", {"IL3RA", "IL5RA", "CSF2RA", "CSF2RB"}),
            ("IL-6 receptor family", {"IL6R", "IL6ST", "IL11RA", "CNTFR", "LIFR", "OSMR", "IL31RA"}),
            ("IL-12 receptor family", {"IL12RB1", "IL12RB2", "IL23R", "IL27RA", "CRLF2"}),
            ("Type II cytokine receptors", {"IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "IFNLR1", "IL10RA", "IL10RB", "IL20RA", "IL20RB", "IL22RA1", "IL22RA2"}),
            ("IL-13 receptors", {"IL13RA1", "IL13RA2"}),
        )
        for name, genes in cytokine_groups:
            if gene in genes:
                return name
        return "Hematopoietin and hormone receptors"
    rtk_groups = (
        ("EGFR/ERBB family", {"EGFR", "ERBB2", "ERBB3", "ERBB4"}),
        ("FGFR family", {"FGFR1", "FGFR2", "FGFR3", "FGFR4"}),
        ("VEGFR family", {"FLT1", "KDR", "FLT4"}),
        ("PDGFR/KIT/CSF1R/FLT3 family", {"PDGFRA", "PDGFRB", "KIT", "CSF1R", "FLT3"}),
        ("Insulin receptor family", {"INSR", "IGF1R", "INSRR"}),
        ("Eph receptors", {gene for gene in _EPH_GENES}),
        ("TAM receptors", {"AXL", "MERTK", "TYRO3"}),
        ("MET family", {"MET", "MST1R"}),
        ("NTRK family", {"NTRK1", "NTRK2", "NTRK3"}),
        ("ALK/LTK family", {"ALK", "LTK"}),
        ("DDR family", {"DDR1", "DDR2"}),
        ("ROR family", {"ROR1", "ROR2"}),
        ("TIE family", {"TEK", "TIE1"}),
    )
    for name, genes in rtk_groups:
        if gene in genes:
            return name
    if gene in {"RET", "ROS1", "RYK", "PTK7", "MUSK", "STYK1"}:
        return f"{gene} receptor family"
    return "Other receptor tyrosine kinases"


_EPH_GENES = {
    *(f"EPHA{value}" for value in (*range(1, 9), 10)),
    *(f"EPHB{value}" for value in range(1, 7)),
}


def _latest_by_gene(rows: Iterable[sqlite3.Row]) -> dict[str, sqlite3.Row]:
    """Role: Compute latest by gene for the global_analysis workflow."""
    output: dict[str, sqlite3.Row] = {}
    for row in rows:
        output.setdefault(str(row["gene"]), row)
    return output


def _connect_source(path: Path) -> sqlite3.Connection:
    """Role: Compute connect source for the global_analysis workflow."""
    connection = sqlite3.connect(f"{path.as_uri()}?mode=ro", uri=True, timeout=2.0)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA query_only = ON")
    return connection


def _has_execution_cache(connection: sqlite3.Connection) -> bool:
    """Role: Compute has execution cache for the global_analysis workflow."""
    return connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='execution_cache'"
    ).fetchone() is not None


def _unique_paths(values: Iterable[str | Path]) -> tuple[Path, ...]:
    """Role: Compute unique paths for the global_analysis workflow."""
    return tuple(dict.fromkeys(Path(value).expanduser().resolve() for value in values))


def _json_object(value: Any) -> dict[str, Any]:
    """Role: Compute json object for the global_analysis workflow."""
    payload = json.loads(str(value)) if isinstance(value, (str, bytes)) else value
    if not isinstance(payload, dict):
        raise TypeError("The global snapshot must be a JSON object")
    return payload


def _records(value: Any) -> list[Mapping[str, Any]]:
    """Role: Compute records for the global_analysis workflow."""
    if not isinstance(value, list):
        return []
    return [item for item in value if isinstance(item, dict)]


def _strings(value: Any) -> tuple[str, ...]:
    """Role: Compute strings for the global_analysis workflow."""
    if not isinstance(value, (list, tuple)):
        return ()
    return tuple(str(item) for item in value if item is not None and str(item))


def _pairs(value: Any) -> tuple[tuple[Any, Any], ...]:
    """Role: Compute pairs for the global_analysis workflow."""
    if not isinstance(value, (list, tuple)):
        return ()
    return tuple(
        (item[0], item[1])
        for item in value
        if isinstance(item, (list, tuple)) and len(item) == 2
    )


def _integer(value: Any) -> int:
    """Role: Compute integer for the global_analysis workflow."""
    try:
        parsed = int(value or 0)
    except (TypeError, ValueError):
        return 0
    return max(0, parsed)


def _hash_object(value: Any) -> str:
    """Role: Compute a stable hash for object."""
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _canonical_json(value: Any) -> str:
    """Role: Compute canonical json for the global_analysis workflow."""
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _utc_now() -> str:
    """Role: Return the current UTC timestamp in canonical form."""
    return datetime.now(timezone.utc).isoformat()


def _write_dataclass_csv(path: Path, records: Iterable[Any]) -> None:
    """Role: Compute write dataclass csv for the global_analysis workflow."""
    rows = [asdict(record) for record in records]
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _safe_csv(value) for key, value in row.items()})


def _safe_csv(value: Any) -> Any:
    """Role: Sanitize csv for safe serialization."""
    if value is None or isinstance(value, (int, float)):
        return value
    text = str(value)
    if text.lstrip(" \t\r\n").startswith(("=", "+", "-", "@")):
        return "'" + text
    return text
