"""Cross-gene analytical projections built from immutable 2.1 snapshots."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .catalog import ExecutionSnapshot
from .disease_categories import snapshot_disease_partition


ACTIVE_REGULATORY_STATUSES = {
    "regulatory_approved",
    "label_biomarker",
    "conditional_approval",
}
ACTIONABILITY_TIERS = (
    "Variant-specific evidence",
    "Active regulatory evidence",
    "Clinical trial",
    "Mechanism target",
    "No therapy evidence",
)


@dataclass(frozen=True)
class FamilyAssignment:
    """Assign one gene to a curated receptor group and protein family."""

    gene: str
    family_group: str
    protein_family: str


@dataclass(frozen=True)
class GlobalGeneSummary:
    """Compact analytical projection of the newest compatible gene snapshots."""

    gene: str
    family_group: str
    protein_family: str
    analysis_execution_id: int
    therapy_execution_id: int | None
    rare_variant_count: int
    cancer_observation_count: int
    shared_variant_count: int
    confirmed_shared_count: int
    candidate_shared_count: int
    max_variant_patient_count: int
    rare_disease_count: int
    hereditary_cancer_disease_count: int
    other_inherited_disease_count: int
    cancer_project_count: int
    variant_therapy_count: int
    active_regulatory_count: int
    clinical_trial_count: int
    mechanism_target_count: int
    actionability_tier: str
    source_coverage_count: int

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-safe row for tables and Plotly figures."""

        return asdict(self)


@dataclass(frozen=True)
class GlobalFamilySummary:
    """Aggregate gene evidence within one curated protein family."""

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
    hereditary_cancer_disease_count: int
    other_inherited_disease_count: int
    cancer_project_count: int

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-safe family row for dashboard components."""

        return asdict(self)


@dataclass(frozen=True)
class GlobalOverview:
    """Global dashboard payload derived without querying external sources."""

    analysis_snapshot_count: int
    active_gene_count: int
    matched_therapy_count: int
    unmatched_therapy_gene_count: int
    shared_variant_count: int
    confirmed_shared_count: int
    actionable_gene_count: int
    source_database_count: int
    gene_summaries: tuple[GlobalGeneSummary, ...]
    family_summaries: tuple[GlobalFamilySummary, ...]
    disease_family_counts: tuple[tuple[str, str, int], ...]
    hereditary_cancer_family_counts: tuple[tuple[str, str, int], ...]
    other_inherited_disease_family_counts: tuple[tuple[str, str, int], ...]
    project_family_counts: tuple[tuple[str, str, int], ...]
    caveat: str = (
        "Evidence layers are displayed separately: variant-specific evidence, "
        "regulatory status, clinical trials and mechanism targets are not equivalent."
    )

    def to_dict(self) -> dict[str, Any]:
        """Return the complete global payload for download or caching."""

        return asdict(self)


def build_global_overview(
    analysis_snapshots: Sequence[ExecutionSnapshot],
    therapy_snapshots: Sequence[ExecutionSnapshot],
    family_catalog: Mapping[str, FamilyAssignment],
) -> GlobalOverview:
    """Project the newest analysis per gene and its compatible therapy snapshot."""

    latest_analyses = _latest_by_gene(analysis_snapshots)
    therapies_by_gene: dict[str, list[ExecutionSnapshot]] = defaultdict(list)
    for snapshot in therapy_snapshots:
        therapies_by_gene[snapshot.gene.strip().upper()].append(snapshot)

    genes: list[GlobalGeneSummary] = []
    disease_family: Counter[tuple[str, str]] = Counter()
    hereditary_cancer_family: Counter[tuple[str, str]] = Counter()
    other_inherited_family: Counter[tuple[str, str]] = Counter()
    project_family: Counter[tuple[str, str]] = Counter()
    matched_therapy_count = 0
    unmatched_therapy_gene_count = 0
    databases: set[Path] = set()

    for gene, analysis_snapshot in sorted(latest_analyses.items()):
        databases.add(analysis_snapshot.database_path)
        result = analysis_snapshot.result
        assignment = family_catalog.get(
            gene,
            FamilyAssignment(gene, "Other analyzed proteins", "Unclassified"),
        )
        analysis_signature = _analysis_signature(result)
        matching_therapy = next(
            (
                snapshot
                for snapshot in therapies_by_gene.get(gene, ())
                if _therapy_signature(snapshot.parameters) == analysis_signature
            ),
            None,
        )
        if matching_therapy is not None:
            matched_therapy_count += 1
            databases.add(matching_therapy.database_path)
            therapy = matching_therapy.result
        else:
            therapy = {}
            if therapies_by_gene.get(gene):
                unmatched_therapy_gene_count += 1

        shared = _records(result.get("shared_mutations"))
        confirmed = sum(
            str(record.get("match_level")) == "confirmed_genomic_and_protein"
            for record in shared
        )
        max_patients = max(
            (_integer(record.get("cancer_patient_count")) for record in shared),
            default=0,
        )
        diseases = {
            disease
            for record in shared
            for disease in _strings(record.get("diseases"))
        }
        hereditary_cancer_diseases = {
            disease
            for record in shared
            for disease in snapshot_disease_partition(record)[0]
        }
        other_inherited_diseases = {
            disease
            for record in shared
            for disease in snapshot_disease_partition(record)[1]
        }
        projects = {
            project
            for record in shared
            for project in _strings(record.get("tcga_projects"))
        }
        for record in shared:
            for disease in _strings(record.get("diseases")):
                disease_family[(disease, assignment.protein_family)] += 1
            hereditary_cancer, other_inherited = snapshot_disease_partition(record)
            for disease in hereditary_cancer:
                hereditary_cancer_family[(disease, assignment.protein_family)] += 1
            for disease in other_inherited:
                other_inherited_family[(disease, assignment.protein_family)] += 1
            for project in _strings(record.get("tcga_projects")):
                project_family[(project, assignment.protein_family)] += 1

        variant_therapy = len(_records(therapy.get("variant_therapy_evidence")))
        regulatory = _records(therapy.get("regulatory_therapy"))
        discovery = _records(therapy.get("trial_and_target_discovery"))
        active_regulatory = sum(
            str(record.get("therapy_status")) in ACTIVE_REGULATORY_STATUSES
            for record in regulatory
        )
        trials = sum(
            str(record.get("therapy_status")) == "clinical_trial"
            for record in discovery
        )
        mechanisms = sum(
            str(record.get("therapy_status")) == "mechanism_target"
            for record in discovery
        )
        source_coverage = len(
            {
                str(record.get("source"))
                for record in _records(result.get("source_statuses"))
                if str(record.get("status")) in {"ok", "partial"}
            }
        )
        genes.append(
            GlobalGeneSummary(
                gene=gene,
                family_group=assignment.family_group,
                protein_family=assignment.protein_family,
                analysis_execution_id=analysis_snapshot.execution_id,
                therapy_execution_id=(
                    matching_therapy.execution_id if matching_therapy else None
                ),
                rare_variant_count=_integer(result.get("rare_variant_count")),
                cancer_observation_count=_integer(
                    result.get("cancer_observation_count")
                ),
                shared_variant_count=len(shared),
                confirmed_shared_count=confirmed,
                candidate_shared_count=len(shared) - confirmed,
                max_variant_patient_count=max_patients,
                rare_disease_count=len(diseases),
                hereditary_cancer_disease_count=len(hereditary_cancer_diseases),
                other_inherited_disease_count=len(other_inherited_diseases),
                cancer_project_count=len(projects),
                variant_therapy_count=variant_therapy,
                active_regulatory_count=active_regulatory,
                clinical_trial_count=trials,
                mechanism_target_count=mechanisms,
                actionability_tier=_actionability_tier(
                    variant_therapy,
                    active_regulatory,
                    trials,
                    mechanisms,
                ),
                source_coverage_count=source_coverage,
            )
        )

    family_summaries = _family_summaries(
        genes,
        disease_family,
        hereditary_cancer_family,
        other_inherited_family,
        project_family,
    )
    return GlobalOverview(
        analysis_snapshot_count=len(analysis_snapshots),
        active_gene_count=len(genes),
        matched_therapy_count=matched_therapy_count,
        unmatched_therapy_gene_count=unmatched_therapy_gene_count,
        shared_variant_count=sum(item.shared_variant_count for item in genes),
        confirmed_shared_count=sum(item.confirmed_shared_count for item in genes),
        actionable_gene_count=sum(
            item.actionability_tier != "No therapy evidence" for item in genes
        ),
        source_database_count=len(databases),
        gene_summaries=tuple(genes),
        family_summaries=tuple(family_summaries),
        disease_family_counts=tuple(
            (disease, family, count)
            for (disease, family), count in sorted(disease_family.items())
        ),
        hereditary_cancer_family_counts=tuple(
            (disease, family, count)
            for (disease, family), count in sorted(
                hereditary_cancer_family.items()
            )
        ),
        other_inherited_disease_family_counts=tuple(
            (disease, family, count)
            for (disease, family), count in sorted(other_inherited_family.items())
        ),
        project_family_counts=tuple(
            (project, family, count)
            for (project, family), count in sorted(project_family.items())
        ),
    )


def load_family_catalog(path: str | Path) -> dict[str, FamilyAssignment]:
    """Load curated receptor-family membership from the 2.1 gene lists."""

    root = Path(path)
    groups = (
        ("Receptor tyrosine kinases", "receptor_tyrosine_kinases.txt"),
        ("TGF-beta receptor superfamily", "tgfb_receptors.txt"),
        ("Kinase-linked cytokine receptors", "kinase_linked_cytokine_receptors.txt"),
    )
    catalog: dict[str, FamilyAssignment] = {}
    for group, filename in groups:
        source = root / filename
        if not source.is_file():
            continue
        for value in source.read_text(encoding="utf-8").splitlines():
            gene = value.strip().upper()
            if not gene:
                continue
            if gene in catalog:
                raise ValueError(f"Gene {gene} appears in more than one family group")
            catalog[gene] = FamilyAssignment(
                gene,
                group,
                _protein_subfamily(gene, group),
            )
    return catalog


def convergence_matrix(
    counts: Iterable[tuple[str, str, int]],
    *,
    top_rows: int = 20,
) -> dict[str, Any]:
    """Create a deterministic top-N entity-by-family matrix."""

    if top_rows < 1:
        raise ValueError("top_rows must be positive")
    triples = tuple(counts)
    row_totals: Counter[str] = Counter()
    for row, _family, count in triples:
        row_totals[row] += count
    rows = [
        row
        for row, _count in sorted(
            row_totals.items(), key=lambda item: (-item[1], item[0])
        )[:top_rows]
    ]
    families = sorted(
        {
            family
            for row, family, count in triples
            if row in rows and count > 0
        }
    )
    lookup = {(row, family): count for row, family, count in triples}
    return {
        "rows": rows,
        "families": families,
        "values": [
            [lookup.get((row, family), 0) for family in families] for row in rows
        ],
    }


def _family_summaries(
    genes: Sequence[GlobalGeneSummary],
    disease_family: Counter[tuple[str, str]],
    hereditary_cancer_family: Counter[tuple[str, str]],
    other_inherited_family: Counter[tuple[str, str]],
    project_family: Counter[tuple[str, str]],
) -> list[GlobalFamilySummary]:
    """Aggregate gene summaries without conflating therapy evidence layers."""

    grouped: dict[tuple[str, str], list[GlobalGeneSummary]] = defaultdict(list)
    for gene in genes:
        grouped[(gene.family_group, gene.protein_family)].append(gene)
    output = []
    for (group, family), members in grouped.items():
        output.append(
            GlobalFamilySummary(
                family_group=group,
                protein_family=family,
                gene_count=len(members),
                genes_with_shared_variants=sum(
                    member.shared_variant_count > 0 for member in members
                ),
                shared_variant_count=sum(
                    member.shared_variant_count for member in members
                ),
                confirmed_shared_count=sum(
                    member.confirmed_shared_count for member in members
                ),
                genes_with_variant_therapy=sum(
                    member.variant_therapy_count > 0 for member in members
                ),
                genes_with_active_regulatory=sum(
                    member.active_regulatory_count > 0 for member in members
                ),
                genes_with_clinical_trials=sum(
                    member.clinical_trial_count > 0 for member in members
                ),
                genes_with_mechanism_target=sum(
                    member.mechanism_target_count > 0 for member in members
                ),
                rare_disease_count=sum(
                    1
                    for disease, item_family in disease_family
                    if item_family == family
                ),
                hereditary_cancer_disease_count=sum(
                    1
                    for disease, item_family in hereditary_cancer_family
                    if item_family == family
                ),
                other_inherited_disease_count=sum(
                    1
                    for disease, item_family in other_inherited_family
                    if item_family == family
                ),
                cancer_project_count=sum(
                    1
                    for project, item_family in project_family
                    if item_family == family
                ),
            )
        )
    return sorted(
        output,
        key=lambda item: (
            -item.shared_variant_count,
            -item.genes_with_variant_therapy,
            item.protein_family,
        ),
    )


def _latest_by_gene(
    snapshots: Sequence[ExecutionSnapshot],
) -> dict[str, ExecutionSnapshot]:
    """Select the first newest snapshot for each normalized gene."""

    output: dict[str, ExecutionSnapshot] = {}
    for snapshot in snapshots:
        output.setdefault(snapshot.gene.strip().upper(), snapshot)
    return output


def _analysis_signature(result: Mapping[str, Any]) -> Mapping[str, Any]:
    """Reconstruct the stable molecular signature stored with therapy queries."""

    return {
        "gene": str(result.get("gene") or "").upper(),
        "shared_mutations": [
            {
                "protein_change": record.get("protein_change"),
                "diseases": list(_strings(record.get("diseases"))),
                "orphanet_codes": list(_strings(record.get("orphanet_codes"))),
                "match_level": record.get("match_level"),
            }
            for record in _records(result.get("shared_mutations"))
        ],
    }


def _therapy_signature(parameters: Mapping[str, Any]) -> Mapping[str, Any] | None:
    """Extract the molecular signature captured by a therapy execution."""

    arguments = parameters.get("arguments")
    if isinstance(arguments, dict) and isinstance(
        arguments.get("analysis_signature"), dict
    ):
        return arguments["analysis_signature"]
    signature = parameters.get("analysis_signature")
    return signature if isinstance(signature, dict) else None


def _actionability_tier(
    variant: int,
    regulatory: int,
    trials: int,
    mechanisms: int,
) -> str:
    """Classify a gene by its highest available, non-equivalent evidence layer."""

    if variant:
        return "Variant-specific evidence"
    if regulatory:
        return "Active regulatory evidence"
    if trials:
        return "Clinical trial"
    if mechanisms:
        return "Mechanism target"
    return "No therapy evidence"


def _protein_subfamily(gene: str, group: str) -> str:
    """Map curated receptor symbols to interpretable protein subfamilies."""

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
        ("Eph receptors", _EPH_GENES),
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


def _records(value: object) -> list[Mapping[str, Any]]:
    """Return mapping records from a JSON-like sequence."""

    if not isinstance(value, (list, tuple)):
        return []
    return [record for record in value if isinstance(record, dict)]


def _strings(value: object) -> tuple[str, ...]:
    """Return non-empty strings from a JSON-like sequence."""

    if not isinstance(value, (list, tuple)):
        return ()
    return tuple(str(item) for item in value if str(item))


def _integer(value: object) -> int:
    """Convert a JSON scalar into a safe non-negative integer."""

    try:
        return max(0, int(value))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return 0


def canonical_signature(value: Mapping[str, Any]) -> str:
    """Expose deterministic signature serialization for diagnostics and tests."""

    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
