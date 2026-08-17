"""Optional Matplotlib visualizations for a gene analysis."""

from __future__ import annotations

import sys
import tempfile
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Any

from .models import AnalysisResult
from .plot_data import (
    heatmap_data,
    lollipop_records,
    mutation_category,
    recurrence_records,
    venn_counts,
)
from .http import ApiError
from .protein_annotations import (
    ProteinAnnotation,
    ProteinFeature,
    fetch_protein_annotation,
    inferred_annotation,
)


COLORS = {
    "rare": "#3B82F6",
    "oncology": "#EF4444",
    "shared": "#8B5CF6",
}
MARKERS = {"rare": "o", "oncology": "o", "shared": "o"}


class PlotDependencyError(RuntimeError):
    """Raised when the optional Matplotlib dependency is unavailable."""


@lru_cache(maxsize=1)
def _matplotlib() -> tuple[Any, Any, Any]:
    """Role: Load the optional Matplotlib plotting dependencies."""
    try:
        import matplotlib

        if "matplotlib.pyplot" not in sys.modules:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.patches import Circle
    except ImportError as exc:
        raise PlotDependencyError(
            "Los graficos requieren Matplotlib. Instale con: "
            "python -m pip install -e '.[plots]'"
        ) from exc
    return plt, Line2D, Circle


@lru_cache(maxsize=1)
def _geneview() -> tuple[Any, Any, Any]:
    """Role: Load the optional Geneview plotting dependencies."""
    try:
        import pandas as pd
        from geneview.genometracks import GenomicInterval, lolliplot
    except ImportError as exc:
        raise PlotDependencyError(
            "El lollipop requiere Geneview. Instale con: "
            "python -m pip install -e '.[plots]'"
        ) from exc
    return pd, lolliplot, GenomicInterval


def create_lollipop_plot(
    result: AnalysisResult,
    *,
    max_mutations: int = 50,
    protein_annotation: ProteinAnnotation | None = None,
) -> Any:
    """Role: Create the lollipop plot visualization or result object."""
    _validate_maximum(max_mutations, "max_mutations")
    annotation = _resolve_protein_annotation(result, protein_annotation)
    plt, _line2d, _circle = _matplotlib()
    fig, ax = plt.subplots(figsize=(13, 6), layout="constrained")
    try:
        _draw_lollipop(
            ax,
            result,
            annotation,
            max_mutations=max_mutations,
        )
    except Exception:
        plt.close(fig)
        raise
    return fig


def create_venn_plot(result: AnalysisResult) -> Any:
    """Role: Create the venn plot visualization or result object."""
    plt, _line2d, _circle = _matplotlib()
    fig, ax = plt.subplots(figsize=(8, 6), layout="constrained")
    try:
        _draw_venn(ax, result)
    except Exception:
        plt.close(fig)
        raise
    return fig


def create_recurrence_plot(result: AnalysisResult) -> Any:
    """Role: Create the recurrence plot visualization or result object."""
    plt, _line2d, _circle = _matplotlib()
    fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")
    try:
        _draw_recurrence(ax, result)
    except Exception:
        plt.close(fig)
        raise
    return fig


def create_heatmap_plot(
    result: AnalysisResult, *, max_mutations: int = 30
) -> Any:
    """Role: Create the heatmap plot visualization or result object."""
    _validate_maximum(max_mutations, "max_mutations")
    plt, _line2d, _circle = _matplotlib()
    fig, ax = plt.subplots(figsize=(12, 8), layout="constrained")
    try:
        image = _draw_heatmap(ax, result, max_mutations=max_mutations)
        if image is not None:
            fig.colorbar(image, ax=ax, label="Unique patients")
    except Exception:
        plt.close(fig)
        raise
    return fig


def create_dashboard(
    result: AnalysisResult,
    *,
    max_lollipop_mutations: int = 50,
    max_heatmap_mutations: int = 30,
    protein_annotation: ProteinAnnotation | None = None,
) -> Any:
    """Role: Create the dashboard visualization or result object."""
    _validate_maximum(max_lollipop_mutations, "max_lollipop_mutations")
    _validate_maximum(max_heatmap_mutations, "max_heatmap_mutations")
    annotation = _resolve_protein_annotation(result, protein_annotation)
    plt, _line2d, _circle = _matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(18, 13), layout="constrained")
    try:
        _draw_lollipop(
            axes[0, 0],
            result,
            annotation,
            max_mutations=max_lollipop_mutations,
        )
        _draw_venn(axes[0, 1], result)
        _draw_recurrence(axes[1, 0], result)
        heatmap = _draw_heatmap(
            axes[1, 1],
            result,
            max_mutations=max_heatmap_mutations,
        )
        if heatmap is not None:
            fig.colorbar(
                heatmap,
                ax=axes[1, 1],
                shrink=0.8,
                label="Unique patients",
            )
        fig.suptitle(
            f"{result.gene}: rare diseases and oncology",
            fontsize=18,
            fontweight="bold",
        )
    except Exception:
        plt.close(fig)
        raise
    return fig


def generate_gene_plots(
    result: AnalysisResult,
    output_dir: str | Path,
    *,
    formats: tuple[str, ...] = ("png", "svg"),
    max_lollipop_mutations: int = 50,
    max_heatmap_mutations: int = 30,
    dpi: int = 200,
    protein_annotation: ProteinAnnotation | None = None,
) -> dict[str, tuple[Path, ...]]:
    """Create the four plots plus a dashboard and return all written paths."""

    normalized_formats = tuple(dict.fromkeys(value.lower() for value in formats))
    unsupported = set(normalized_formats) - {"png", "svg"}
    if not normalized_formats or unsupported:
        raise ValueError("formats debe contener solamente 'png' y/o 'svg'")
    if dpi < 72:
        raise ValueError("dpi must be at least 72")
    _validate_maximum(max_lollipop_mutations, "max_lollipop_mutations")
    _validate_maximum(max_heatmap_mutations, "max_heatmap_mutations")
    annotation = _resolve_protein_annotation(result, protein_annotation)

    creators = {
        "lollipop": lambda: create_lollipop_plot(
            result,
            max_mutations=max_lollipop_mutations,
            protein_annotation=annotation,
        ),
        "venn": lambda: create_venn_plot(result),
        "recurrence": lambda: create_recurrence_plot(result),
        "heatmap": lambda: create_heatmap_plot(
            result, max_mutations=max_heatmap_mutations
        ),
        "dashboard": lambda: create_dashboard(
            result,
            max_lollipop_mutations=max_lollipop_mutations,
            max_heatmap_mutations=max_heatmap_mutations,
            protein_annotation=annotation,
        ),
    }
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    plt, _line2d, _circle = _matplotlib()
    written: dict[str, tuple[Path, ...]] = {}
    with tempfile.TemporaryDirectory(prefix=".rare-oncology-", dir=destination) as tmp:
        temporary = Path(tmp)
        staged: dict[str, list[tuple[Path, Path]]] = {}
        for name, creator in creators.items():
            figure = creator()
            paths: list[tuple[Path, Path]] = []
            try:
                for image_format in normalized_formats:
                    filename = f"{result.gene.lower()}_{name}.{image_format}"
                    temporary_path = temporary / filename
                    final_path = destination / filename
                    figure.savefig(
                        temporary_path,
                        dpi=dpi,
                        bbox_inches="tight",
                        metadata={
                            "Creator": "rare-oncology 2.0",
                            "Title": f"{result.gene} {name}",
                        },
                    )
                    paths.append((temporary_path, final_path))
            finally:
                plt.close(figure)
            staged[name] = paths
        for name, paths in staged.items():
            for temporary_path, final_path in paths:
                temporary_path.replace(final_path)
            written[name] = tuple(final_path for _temporary, final_path in paths)
    return written


def _draw_lollipop(
    ax: Any,
    result: AnalysisResult,
    annotation: ProteinAnnotation,
    *,
    max_mutations: int,
) -> None:
    """Role: Draw lollipop on the supplied plotting axis."""
    canonical_limit = annotation.length if annotation.accession else None
    records = lollipop_records(
        result,
        max_mutations=max_mutations,
        max_position=canonical_limit,
    )
    total_positioned = sum(
        record.amino_acid_position is not None
        and (
            canonical_limit is None
            or record.amino_acid_position <= canonical_limit
        )
        for record in result.mutation_plot_data
    )
    out_of_range = sum(
        record.amino_acid_position is not None
        and canonical_limit is not None
        and record.amino_acid_position > canonical_limit
        for record in result.mutation_plot_data
    )
    suffix = (
        f" (top {len(records)} de {total_positioned})"
        if total_positioned > len(records)
        else ""
    )
    provenance = (
        f" · UniProt {annotation.accession}"
        if annotation.accession
        else " · inferred length (no UniProt)"
    )
    range_note = (
        f" · {out_of_range} out-of-range change omitted"
        f"{'s' if out_of_range != 1 else ''}"
        if out_of_range
        else ""
    )
    title = f"{result.gene}: variant protein map{suffix}{provenance}{range_note}"
    if not records:
        ax.set_title(title)
        _empty_axis(ax, "No positionable protein changes")
        return
    pd, lolliplot, GenomicInterval = _geneview()
    maximum_position = max(record.amino_acid_position or 1 for record in records)
    protein_length = (
        annotation.length
        if annotation.accession
        else max(annotation.length, maximum_position)
    )
    variants = pd.DataFrame(
        [
            {
                "chrom": result.gene,
                "start": record.amino_acid_position,
                "end": (record.amino_acid_position or 1) + 1,
                "score": max(record.cancer_patient_count, 1),
                "label": (
                    f"p.{record.protein_change}"
                    if len(records) <= 12
                    or record.in_rare_disease
                    else ""
                ),
                "fill": COLORS[mutation_category(record)],
                "color": COLORS[mutation_category(record)],
                "border": "#FFFFFF",
                "shape": "circle",
                "label_rotation": 0,
                "cex": 1.0,
            }
            for record in records
        ]
    )
    domains = _geneview_domains(pd, result.gene, annotation.domains)
    lolliplot(
        variants,
        features=domains,
        type="circle",
        region=GenomicInterval(result.gene, 1, protein_length + 1),
        ax=ax,
        show_yaxis=False,
        label_on_feature=True,
        lollipop_style_switch_limit=0,
        node_size=0.014,
        cex=1.0,
        jitter="label",
    )
    _spread_geneview_labels(ax, protein_length)
    ax.set_ylim(-0.34, 1.12)
    _draw_secondary_structure(ax, annotation.secondary_structure, protein_length)
    _add_track_labels(ax)
    _add_source_legend(ax)
    ax.set_xlabel("Residue")
    ax.set_xlim(1, protein_length + 1)
    residue_ticks = _residue_ticks(protein_length)
    ax.set_xticks(residue_ticks, labels=[str(value) for value in residue_ticks])
    ax.tick_params(axis="x", labelsize=8)
    ax.text(
        0.99,
        0.01,
        "α  Helix     →  β-strand",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color="#475569",
    )
    ax.set_title(title, pad=34)


def _draw_venn(ax: Any, result: AnalysisResult) -> None:
    """Role: Draw venn on the supplied plotting axis."""
    _plt, _line2d, Circle = _matplotlib()
    counts = venn_counts(result)
    left = Circle(
        (0.42, 0.5),
        0.29,
        facecolor=COLORS["rare"],
        edgecolor="#1E3A8A",
        alpha=0.4,
    )
    right = Circle(
        (0.62, 0.5),
        0.29,
        facecolor=COLORS["oncology"],
        edgecolor="#991B1B",
        alpha=0.4,
    )
    ax.add_patch(left)
    ax.add_patch(right)
    ax.text(0.29, 0.5, str(counts.rare_only), ha="center", va="center", fontsize=18)
    ax.text(0.52, 0.5, str(counts.shared), ha="center", va="center", fontsize=18)
    ax.text(0.75, 0.5, str(counts.oncology_only), ha="center", va="center", fontsize=18)
    rare_total = counts.rare_only + counts.shared
    oncology_total = counts.oncology_only + counts.shared
    ax.text(
        0.28,
        0.17,
        f"ClinVar–Orphanet (n={rare_total})",
        ha="center",
        fontsize=11,
    )
    ax.text(
        0.76,
        0.17,
        f"GDC/TCGA (n={oncology_total})",
        ha="center",
        fontsize=11,
    )
    ax.text(
        0.52,
        0.09,
        "Schematic diagram: areas are not proportional",
        ha="center",
        fontsize=8,
        color="#475569",
    )
    ax.set_title(f"{result.gene}: unique protein variants")
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0.08, 0.9)
    ax.set_aspect("equal")
    ax.axis("off")


def _draw_recurrence(ax: Any, result: AnalysisResult) -> None:
    """Role: Draw recurrence on the supplied plotting axis."""
    mutations = recurrence_records(result)
    ax.set_title(f"{result.gene}: shared-mutation recurrence")
    if not mutations:
        _empty_axis(ax, "No shared mutations under the current filters")
        return
    labels = [mutation.protein_change for mutation in mutations]
    values = [mutation.cancer_patient_count for mutation in mutations]
    colors = [
        "#16A34A"
        if mutation.match_level == "confirmed_genomic_and_protein"
        else COLORS["shared"]
        for mutation in mutations
    ]
    bars = ax.bar(labels, values, color=colors)
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.set_xlabel("Cambio proteico")
    ax.set_ylabel("Unique patients")
    ax.tick_params(axis="x", rotation=45)
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    _add_evidence_legend(ax)


def _draw_heatmap(
    ax: Any,
    result: AnalysisResult,
    *,
    max_mutations: int,
) -> Any | None:
    """Role: Draw heatmap on the supplied plotting axis."""
    data = heatmap_data(result, max_mutations=max_mutations)
    total_oncology = sum(
        record.in_oncology for record in result.mutation_plot_data
    )
    suffix = (
        f" (top {len(data.mutations)} de {total_oncology})"
        if total_oncology > len(data.mutations)
        else ""
    )
    ax.set_title(f"{result.gene}: patients by mutation and project{suffix}")
    if not data.mutations or not data.projects:
        _empty_axis(ax, "No oncology observations are available for the heatmap")
        return None
    image = ax.imshow(
        data.values,
        cmap="Blues",
        interpolation="nearest",
        aspect="auto",
        vmin=0,
    )
    ax.set_xticks(range(len(data.projects)), labels=data.projects, rotation=45, ha="right")
    ax.set_yticks(range(len(data.mutations)), labels=data.mutations)
    ax.set_xlabel("Proyecto GDC")
    ax.set_ylabel("Cambio proteico")
    if len(data.mutations) * len(data.projects) <= 120:
        maximum = max(max(row) for row in data.values) or 1
        for row_index, row in enumerate(data.values):
            for column_index, value in enumerate(row):
                ax.text(
                    column_index,
                    row_index,
                    str(value),
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if value > maximum / 2 else "#0F172A",
                )
    return image


def _add_source_legend(ax: Any) -> None:
    """Role: Add source legend annotations to the current plot."""
    _plt, Line2D, _circle = _matplotlib()
    handles = [
        Line2D(
            [0],
            [0],
            marker=MARKERS[category],
            color="none",
            markerfacecolor=COLORS[category],
            markeredgecolor="white",
            markersize=9,
            label=label,
        )
        for category, label in (
            ("rare", "Rare disease"),
            ("oncology", "Oncology"),
            ("shared", "Shared"),
        )
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8)


def _add_evidence_legend(ax: Any) -> None:
    """Role: Add evidence legend annotations to the current plot."""
    _plt, Line2D, _circle = _matplotlib()
    handles = [
        Line2D([0], [0], color="#16A34A", linewidth=8, label="Genomically confirmed"),
        Line2D([0], [0], color=COLORS["shared"], linewidth=8, label="Candidate/mixed"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8)


def _empty_axis(ax: Any, message: str) -> None:
    """Role: Render an empty-state axis placeholder."""
    ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def _validate_maximum(value: int, name: str) -> None:
    """Role: Validate maximum constraints before execution."""
    if value < 1:
        raise ValueError(f"{name} debe ser positivo")


def _resolve_protein_annotation(
    result: AnalysisResult,
    provided: ProteinAnnotation | None,
) -> ProteinAnnotation:
    """Role: Resolve protein annotation from cached, supplied, or inferred data."""
    if provided is not None:
        return provided
    positions = [
        record.amino_acid_position
        for record in result.mutation_plot_data
        if record.amino_acid_position is not None
    ]
    lookup_failed = False
    try:
        annotation = fetch_protein_annotation(result.gene)
    except ApiError as exc:
        lookup_failed = True
        warnings.warn(
            f"No fue posible consultar UniProt para {result.gene}: {exc}",
            RuntimeWarning,
            stacklevel=2,
        )
        annotation = None
    if annotation is None and not lookup_failed:
        warnings.warn(
            f"No hay una entrada humana revisada de UniProt para {result.gene}; "
            "an inferred length without domains or secondary structure will be used.",
            RuntimeWarning,
            stacklevel=2,
        )
    return annotation or inferred_annotation(result.gene, positions)


def _geneview_domains(pd: Any, gene: str, domains: tuple[ProteinFeature, ...]) -> Any:
    """Role: Load the optional Geneview plotting dependencies."""
    if not domains:
        return None
    palette = ("#BFDBFE", "#C4B5FD", "#A7F3D0", "#FDE68A", "#FBCFE8")
    return pd.DataFrame(
        [
            {
                "chrom": gene,
                "start": domain.start,
                "end": domain.end + 1,
                "name": domain.label,
                "fill": palette[index % len(palette)],
                "color": "#475569",
                "height": 0.065,
                "feature": "domain",
            }
            for index, domain in enumerate(domains)
        ]
    )


def _draw_secondary_structure(
    ax: Any,
    features: tuple[ProteinFeature, ...],
    protein_length: int,
) -> None:
    """Role: Draw secondary structure on the supplied plotting axis."""
    y = -0.16
    ax.hlines(y, 1, protein_length + 1, color="#94A3B8", linewidth=0.8, zorder=1)
    if not features:
        ax.text(
            protein_length / 2,
            y,
            "No secondary-structure annotation is available in UniProt",
            ha="center",
            va="center",
            fontsize=7,
            color="#64748B",
        )
        return
    for feature in features:
        if feature.kind == "helix":
            ax.plot(
                [feature.start, feature.end + 1],
                [y, y],
                color="#F59E0B",
                linewidth=7,
                solid_capstyle="round",
                zorder=3,
            )
            if feature.end - feature.start >= 5:
                alpha_count = min(6, max(1, (feature.end - feature.start) // 4))
                ax.text(
                    (feature.start + feature.end) / 2,
                    y,
                    "α" * alpha_count,
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="white",
                    zorder=4,
                )
        elif feature.kind == "strand":
            ax.annotate(
                "",
                xy=(feature.end + 1, y),
                xytext=(feature.start, y),
                arrowprops={
                    "arrowstyle": "-|>",
                    "color": "#0EA5E9",
                    "linewidth": 3,
                    "mutation_scale": 8,
                },
                zorder=3,
            )


def _add_track_labels(ax: Any) -> None:
    """Role: Add track labels annotations to the current plot."""
    transform = ax.get_yaxis_transform()
    for y, label in ((0.48, "Variants"), (0.08, "Domains"), (-0.16, "Sec. structure")):
        ax.text(
            -0.015,
            y,
            label,
            transform=transform,
            ha="right",
            va="center",
            fontsize=9,
            fontweight="bold",
            clip_on=False,
        )


def _spread_geneview_labels(ax: Any, protein_length: int) -> None:
    """Role: Compute spread geneview labels for the plots workflow."""
    labels = sorted(
        (text for text in ax.texts if text.get_text().startswith("p.")),
        key=lambda text: text.get_position()[0],
    )
    if len(labels) < 2:
        return
    levels = (0.86, 0.94, 1.02, 1.10)
    last_x = [float("-inf")] * len(levels)
    minimum_gap = max(12.0, protein_length * 0.035)
    for label in labels:
        x, _y = label.get_position()
        level_index = next(
            (
                index
                for index, previous_x in enumerate(last_x)
                if x - previous_x >= minimum_gap
            ),
            min(range(len(levels)), key=last_x.__getitem__),
        )
        label.set_position((x, levels[level_index]))
        last_x[level_index] = x


def _residue_ticks(protein_length: int, intervals: int = 5) -> tuple[int, ...]:
    """Role: Compute residue ticks for the plots workflow."""
    if protein_length <= 1:
        return (1,)
    return tuple(
        sorted(
            {
                1,
                protein_length,
                *(
                    round(1 + (protein_length - 1) * index / intervals)
                    for index in range(1, intervals)
                ),
            }
        )
    )
