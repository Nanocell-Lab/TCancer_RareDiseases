"""Pure data transforms used by the optional plotting layer."""

from __future__ import annotations

from dataclasses import dataclass

from .models import AnalysisResult, MutationPlotRecord, SharedMutation


@dataclass(frozen=True)
class VennCounts:
    rare_only: int
    shared: int
    oncology_only: int


@dataclass(frozen=True)
class HeatmapData:
    mutations: tuple[str, ...]
    projects: tuple[str, ...]
    values: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class RareDiseaseHeatmapData:
    mutations: tuple[str, ...]
    diseases: tuple[str, ...]
    values: tuple[tuple[int, ...], ...]
    shared_mutations: tuple[bool, ...]


@dataclass(frozen=True)
class LollipopResidue:
    position: int
    protein_changes: tuple[str, ...]
    rare_count: int
    oncology_count: int
    shared_count: int
    cancer_patient_count: int
    confirmed_genomic: bool
    rare_evidence_priority: int
    label_change: str

    @property
    def category(self) -> str:
        """Role: Provide the category operation for the plot_data workflow."""
        if self.shared_count:
            return "shared"
        if self.oncology_count:
            return "oncology"
        return "rare"


def mutation_category(record: MutationPlotRecord) -> str:
    """Role: Provide the mutation category operation for the plot_data workflow."""
    if record.in_rare_disease and record.in_oncology:
        return "shared"
    if record.in_rare_disease:
        return "rare"
    return "oncology"


def venn_counts(result: AnalysisResult) -> VennCounts:
    """Role: Provide the venn counts operation for the plot_data workflow."""
    rare_only = 0
    shared = 0
    oncology_only = 0
    for record in result.mutation_plot_data:
        category = mutation_category(record)
        if category == "rare":
            rare_only += 1
        elif category == "shared":
            shared += 1
        else:
            oncology_only += 1
    return VennCounts(rare_only, shared, oncology_only)


def lollipop_records(
    result: AnalysisResult,
    *,
    max_mutations: int = 50,
    max_position: int | None = None,
) -> tuple[MutationPlotRecord, ...]:
    """Role: Provide the lollipop records operation for the plot_data workflow."""
    if max_mutations < 1:
        raise ValueError("max_mutations must be positive")
    positioned = [
        record
        for record in result.mutation_plot_data
        if record.amino_acid_position is not None
        and (max_position is None or record.amino_acid_position <= max_position)
    ]
    positioned.sort(
        key=lambda record: (
            {"shared": 0, "rare": 1, "oncology": 2}[mutation_category(record)],
            -record.cancer_patient_count,
            record.amino_acid_position or 0,
            record.protein_change,
        )
    )
    selected = positioned[:max_mutations]
    return tuple(
        sorted(
            selected,
            key=lambda record: (
                record.amino_acid_position or 0,
                record.protein_change,
            ),
        )
    )


def lollipop_residues(
    result: AnalysisResult,
    *,
    max_residues: int = 80,
    max_position: int | None = None,
    include_synonymous: bool = False,
) -> tuple[LollipopResidue, ...]:
    """Return an exact-position, category-balanced protein overview.

    Multiple changes at one amino-acid position are aggregated into one node.
    Selection cycles through shared, oncology-only and rare-only residues so a
    large ClinVar set cannot crowd all oncology hotspots out of the plot.
    """

    if max_residues < 1:
        raise ValueError("max_residues must be positive")
    grouped: dict[int, list[MutationPlotRecord]] = {}
    for record in result.mutation_plot_data:
        position = record.amino_acid_position
        if position is None or (max_position is not None and position > max_position):
            continue
        if not include_synonymous and record.protein_change.endswith("="):
            continue
        grouped.setdefault(position, []).append(record)

    residues = tuple(
        _aggregate_residue(position, records)
        for position, records in sorted(grouped.items())
    )
    buckets = {
        category: sorted(
            (residue for residue in residues if residue.category == category),
            key=_lollipop_priority,
        )
        for category in ("shared", "oncology", "rare")
    }
    selected: list[LollipopResidue] = []
    pattern = ("shared", "oncology", "shared", "oncology", "rare")
    indices = {category: 0 for category in buckets}
    while len(selected) < max_residues:
        progressed = False
        for category in pattern:
            index = indices[category]
            if index >= len(buckets[category]):
                continue
            selected.append(buckets[category][index])
            indices[category] += 1
            progressed = True
            if len(selected) == max_residues:
                break
        if not progressed:
            break
    return tuple(sorted(selected, key=lambda residue: residue.position))


def lollipop_label_residues(
    residues: tuple[LollipopResidue, ...], *, max_labels: int = 15
) -> tuple[LollipopResidue, ...]:
    """Choose a small balanced subset of nodes to label."""

    if max_labels < 0:
        raise ValueError("max_labels cannot be negative")
    if max_labels == 0:
        return ()
    buckets = {
        category: sorted(
            (residue for residue in residues if residue.category == category),
            key=_lollipop_priority,
        )
        for category in ("shared", "oncology", "rare")
    }
    pattern = ("shared", "oncology", "shared", "oncology", "rare")
    indices = {category: 0 for category in buckets}
    selected: list[LollipopResidue] = []
    while len(selected) < max_labels:
        progressed = False
        for category in pattern:
            index = indices[category]
            if index >= len(buckets[category]):
                continue
            selected.append(buckets[category][index])
            indices[category] += 1
            progressed = True
            if len(selected) == max_labels:
                break
        if not progressed:
            break
    return tuple(sorted(selected, key=lambda residue: residue.position))


def _aggregate_residue(
    position: int, records: list[MutationPlotRecord]
) -> LollipopResidue:
    """Role: Compute aggregate residue for the plot_data workflow."""
    ordered = sorted(records, key=_record_priority)
    return LollipopResidue(
        position=position,
        protein_changes=tuple(record.protein_change for record in ordered),
        rare_count=sum(mutation_category(record) == "rare" for record in records),
        oncology_count=sum(
            mutation_category(record) == "oncology" for record in records
        ),
        shared_count=sum(mutation_category(record) == "shared" for record in records),
        cancer_patient_count=max(
            (record.cancer_patient_count for record in records), default=0
        ),
        confirmed_genomic=any(
            record.match_level == "confirmed_genomic_and_protein"
            for record in records
        ),
        rare_evidence_priority=min(
            (
                _classification_priority(record.germline_classifications)
                for record in records
                if record.in_rare_disease
            ),
            default=9,
        ),
        label_change=ordered[0].protein_change,
    )


def _record_priority(record: MutationPlotRecord) -> tuple[int, int, int, int, str]:
    """Role: Compute record priority for the plot_data workflow."""
    category = mutation_category(record)
    return (
        0 if record.match_level == "confirmed_genomic_and_protein" else 1,
        {"shared": 0, "oncology": 1, "rare": 2}[category],
        -record.cancer_patient_count,
        _classification_priority(record.germline_classifications),
        record.protein_change,
    )


def _lollipop_priority(residue: LollipopResidue) -> tuple[int, int, int, int, str]:
    """Role: Compute lollipop priority for the plot_data workflow."""
    return (
        0 if residue.confirmed_genomic else 1,
        -residue.cancer_patient_count,
        residue.rare_evidence_priority,
        residue.position,
        residue.label_change,
    )


def _classification_priority(classifications: tuple[str, ...]) -> int:
    """Role: Compute classification priority for the plot_data workflow."""
    text = " | ".join(classifications).casefold()
    if "pathogenic" in text and "conflict" not in text:
        return 0 if "likely" not in text else 1
    if "conflict" in text:
        return 2
    if "uncertain" in text:
        return 3
    if "likely benign" in text:
        return 5
    if "benign" in text:
        return 6
    return 4


def recurrence_records(result: AnalysisResult) -> tuple[SharedMutation, ...]:
    """Role: Provide the recurrence records operation for the plot_data workflow."""
    return tuple(
        sorted(
            result.shared_mutations,
            key=lambda mutation: (
                -mutation.cancer_patient_count,
                mutation.protein_change,
            ),
        )
    )


def heatmap_data(
    result: AnalysisResult, *, max_mutations: int = 30
) -> HeatmapData:
    """Role: Provide the heatmap data operation for the plot_data workflow."""
    if max_mutations < 1:
        raise ValueError("max_mutations must be positive")
    oncology = [
        record for record in result.mutation_plot_data if record.in_oncology
    ]
    oncology.sort(
        key=lambda record: (-record.cancer_patient_count, record.protein_change)
    )
    oncology = oncology[:max_mutations]
    projects = tuple(
        sorted(
            {
                project
                for record in oncology
                for project, _count in record.project_patient_counts
            }
        )
    )
    values = tuple(
        tuple(dict(record.project_patient_counts).get(project, 0) for project in projects)
        for record in oncology
    )
    return HeatmapData(
        mutations=tuple(record.protein_change for record in oncology),
        projects=projects,
        values=values,
    )


def rare_disease_heatmap_data(
    result: AnalysisResult, *, max_mutations: int = 30
) -> RareDiseaseHeatmapData:
    """Build a rare-disease by protein-change association matrix.

    Shared cancer mutations are selected before rare-only mutations. Within
    each group, records with more ClinVar disease associations rank first.
    Older cached results without per-disease counts fall back to one normalized
    association for every disease recorded on the mutation.
    """

    if not 1 <= max_mutations <= 30:
        raise ValueError("max_mutations must be between 1 and 30")
    rare_records = [
        record for record in result.mutation_plot_data if record.in_rare_disease
    ]

    def counts(record: MutationPlotRecord) -> dict[str, int]:
        """Role: Provide the counts operation for the plot_data workflow."""
        if record.disease_variant_counts:
            return {
                disease: count
                for disease, count in record.disease_variant_counts
                if disease and count > 0
            }
        return {disease: 1 for disease in record.diseases if disease}

    def priority(record: MutationPlotRecord) -> tuple[int, int, int, str]:
        """Role: Provide the priority operation for the plot_data workflow."""
        disease_counts = counts(record)
        return (
            -sum(disease_counts.values()),
            -len(disease_counts),
            -record.cancer_patient_count,
            record.protein_change,
        )

    shared = sorted(
        (record for record in rare_records if record.in_oncology), key=priority
    )
    rare_only = sorted(
        (record for record in rare_records if not record.in_oncology), key=priority
    )
    selected = (shared + rare_only)[:max_mutations]
    selected_counts = {
        record.protein_change: counts(record) for record in selected
    }
    diseases = tuple(
        sorted(
            {
                disease
                for disease_counts in selected_counts.values()
                for disease in disease_counts
            },
            key=lambda disease: (
                -sum(
                    disease_counts.get(disease, 0)
                    for disease_counts in selected_counts.values()
                ),
                disease,
            ),
        )
    )
    return RareDiseaseHeatmapData(
        mutations=tuple(record.protein_change for record in selected),
        diseases=diseases,
        values=tuple(
            tuple(
                selected_counts[record.protein_change].get(disease, 0)
                for record in selected
            )
            for disease in diseases
        ),
        shared_mutations=tuple(record.in_oncology for record in selected),
    )
