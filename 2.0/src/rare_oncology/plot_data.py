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
        raise ValueError("max_mutations debe ser positivo")
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
        raise ValueError("max_mutations debe ser positivo")
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
