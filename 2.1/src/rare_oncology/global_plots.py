"""Static global figures generated from the incremental family report."""

from __future__ import annotations

import math
import sys
import tempfile
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

from .global_analysis import GlobalFamilyReport
from .plots import PlotDependencyError


TIER_ORDER = (
    "Active regulatory biomarker",
    "Variant-specific evidence",
    "Clinical investigation",
    "Mechanism/target only",
    "No therapy evidence",
)
TIER_COLORS = {
    "Variant-specific evidence": "#7C3AED",
    "Active regulatory biomarker": "#DC2626",
    "Clinical investigation": "#F59E0B",
    "Mechanism/target only": "#0EA5E9",
    "No therapy evidence": "#CBD5E1",
}
GROUP_COLORS = {
    "Receptor tyrosine kinases": "#2563EB",
    "TGF-beta receptor superfamily": "#7C3AED",
    "Kinase-linked cytokine receptors": "#059669",
    "Other analyzed proteins": "#64748B",
}


@lru_cache(maxsize=1)
def _matplotlib() -> tuple[Any, Any]:
    """Role: Load the optional Matplotlib plotting dependencies."""
    try:
        import matplotlib  # pyright: ignore[reportMissingImports]

        if "matplotlib.pyplot" not in sys.modules:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # pyright: ignore[reportMissingImports]
        from matplotlib.patches import Patch  # pyright: ignore[reportMissingImports]
    except ImportError as exc:
        raise PlotDependencyError(
            "Global plots require Matplotlib. Install it with: "
            "python -m pip install -e '.[plots]'"
        ) from exc
    return plt, Patch


def create_family_actionability_heatmap(report: GlobalFamilyReport) -> Any:
    """Return a normalized family-by-evidence heatmap."""

    plt, _patch = _matplotlib()
    families = list(report.family_summaries[:24])
    fig, ax = plt.subplots(
        figsize=(12, max(5.5, 0.42 * len(families) + 2.2)), layout="constrained"
    )
    if not families:
        _empty_panel(ax, "Protein-family actionability matrix")
        return fig
    columns = (
        "Confirmed\noverlap",
        "Variant-specific\ntherapy",
        "Regulatory\nbiomarker",
        "Clinical\ntrial",
        "Mechanism /\ntarget",
    )
    matrix: list[list[float]] = []
    for family in families:
        confirmed_denominator = max(1, family.shared_variant_count)
        gene_denominator = max(1, family.gene_count)
        matrix.append(
            [
                100 * family.confirmed_shared_count / confirmed_denominator,
                100 * family.genes_with_variant_therapy / gene_denominator,
                100 * family.genes_with_active_regulatory / gene_denominator,
                100 * family.genes_with_clinical_trials / gene_denominator,
                100 * family.genes_with_mechanism_target / gene_denominator,
            ]
        )
    image = ax.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(columns)), columns)
    ax.set_yticks(
        range(len(families)),
        [f"{item.protein_family}  (n={item.gene_count})" for item in families],
    )
    for row, values in enumerate(matrix):
        for column, value in enumerate(values):
            color = "white" if value >= 58 else "#0F172A"
            ax.text(column, row, f"{value:.0f}%", ha="center", va="center", color=color, fontsize=8)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label("Percentage within family")
    ax.set_title("Protein-family actionability matrix", loc="left", fontweight="bold")
    ax.set_xlabel("Evidence layer (not interchangeable)")
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=8)
    return fig


def create_gene_family_landscape(report: GlobalFamilyReport) -> Any:
    """Return a gene landscape of molecular overlap and cancer recurrence."""

    plt, Patch = _matplotlib()
    genes = [item for item in report.gene_summaries if item.shared_variant_count]
    fig, ax = plt.subplots(figsize=(12, 7), layout="constrained")
    if not genes:
        _empty_panel(ax, "Gene-level rare-oncology landscape")
        return fig
    for group in sorted({item.family_group for item in genes}):
        members = [item for item in genes if item.family_group == group]
        ax.scatter(
            [item.confirmed_shared_count for item in members],
            [item.max_variant_patient_count for item in members],
            s=[35 + 16 * math.sqrt(item.shared_variant_count) for item in members],
            alpha=0.72,
            color=GROUP_COLORS.get(group, "#64748B"),
            edgecolor="white",
            linewidth=0.7,
            label=group,
        )
    ranked = sorted(
        genes,
        key=lambda item: (
            item.variant_specific_therapy_count > 0,
            item.max_variant_patient_count,
            item.confirmed_shared_count,
        ),
        reverse=True,
    )[:16]
    for item in ranked:
        ax.annotate(
            item.gene,
            (item.confirmed_shared_count, item.max_variant_patient_count),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xlabel("Confirmed shared rare-oncology variants")
    ax.set_ylabel("Maximum unique cancer patients for one shared variant")
    ax.set_title("Gene-level rare-oncology landscape", loc="left", fontweight="bold")
    ax.grid(alpha=0.2, linewidth=0.7)
    ax.legend(
        handles=[
            Patch(facecolor=color, label=group)
            for group, color in GROUP_COLORS.items()
            if any(item.family_group == group for item in genes)
        ],
        loc="best",
        frameon=False,
        fontsize=8,
    )
    return fig


def create_family_actionability_stack(report: GlobalFamilyReport) -> Any:
    """Return counts of genes by their highest supported actionability tier."""

    plt, _patch = _matplotlib()
    family_tiers: dict[str, Counter[str]] = {}
    for gene in report.gene_summaries:
        family_tiers.setdefault(gene.protein_family, Counter())[gene.actionability_tier] += 1
    ordered = sorted(
        family_tiers,
        key=lambda family: (
            -family_tiers[family]["Active regulatory biomarker"],
            -family_tiers[family]["Variant-specific evidence"],
            family,
        ),
    )[:24]
    fig, ax = plt.subplots(
        figsize=(12, max(5.5, 0.4 * len(ordered) + 2.2)), layout="constrained"
    )
    if not ordered:
        _empty_panel(ax, "Highest actionability tier by family")
        return fig
    left = [0] * len(ordered)
    for tier in TIER_ORDER:
        values = [family_tiers[family][tier] for family in ordered]
        bars = ax.barh(
            range(len(ordered)), values, left=left, color=TIER_COLORS[tier], label=tier
        )
        for bar, value in zip(bars, values):
            if value:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_y() + bar.get_height() / 2,
                    str(value),
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="white" if tier != "No therapy evidence" else "#334155",
                )
        left = [current + value for current, value in zip(left, values)]
    ax.set_yticks(range(len(ordered)), ordered)
    ax.invert_yaxis()
    ax.set_xlabel("Genes (classified by highest available evidence layer)")
    ax.set_title("Highest actionability tier by protein family", loc="left", fontweight="bold")
    ax.legend(ncol=2, frameon=False, fontsize=8, loc="lower right")
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(axis="x", alpha=0.18)
    return fig


def create_rare_disease_family_heatmap(
    report: GlobalFamilyReport, *, top_diseases: int = 20, top_families: int = 16
) -> Any:
    """Return the most abundant rare-disease/family shared-mutation matrix."""

    if top_diseases < 1 or top_families < 1:
        raise ValueError("top_diseases and top_families must be positive")
    plt, _patch = _matplotlib()
    disease_totals: Counter[str] = Counter()
    family_totals: Counter[str] = Counter()
    cells: Counter[tuple[str, str]] = Counter()
    for disease, family, count in report.disease_family_counts:
        cells[(disease, family)] += count
        disease_totals[disease] += count
        family_totals[family] += count
    diseases = [item for item, _count in disease_totals.most_common(top_diseases)]
    families = [item for item, _count in family_totals.most_common(top_families)]
    fig, ax = plt.subplots(
        figsize=(max(10, 0.75 * len(families) + 4), max(6, 0.4 * len(diseases) + 2.5)),
        layout="constrained",
    )
    if not diseases or not families:
        _empty_panel(ax, "Rare-disease convergence across protein families")
        return fig
    matrix = [[cells[(disease, family)] for family in families] for disease in diseases]
    maximum = max(max(row) for row in matrix)
    image = ax.imshow(matrix, cmap="magma", vmin=0, vmax=max(1, maximum), aspect="auto")
    ax.set_xticks(range(len(families)), families, rotation=45, ha="right")
    ax.set_yticks(range(len(diseases)), diseases)
    if len(diseases) * len(families) <= 260:
        for row, values in enumerate(matrix):
            for column, value in enumerate(values):
                if value:
                    ax.text(column, row, str(value), ha="center", va="center", fontsize=7, color="white")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label("Cancer-shared mutation associations")
    ax.set_title(
        f"Top {len(diseases)} rare diseases across protein families",
        loc="left",
        fontweight="bold",
    )
    ax.set_xlabel("Protein family")
    ax.set_ylabel("Rare disease")
    ax.tick_params(labelsize=8)
    return fig


def generate_global_family_plots(
    report: GlobalFamilyReport,
    output_dir: str | Path,
    *,
    formats: Iterable[str] = ("png",),
    dpi: int = 200,
) -> tuple[Path, ...]:
    """Write all implemented global figures and return their paths."""

    normalized_formats = tuple(dict.fromkeys(value.lower().lstrip(".") for value in formats))
    if not normalized_formats or any(value not in {"png", "svg", "pdf"} for value in normalized_formats):
        raise ValueError("formats must contain png, svg, or pdf")
    if dpi < 72:
        raise ValueError("dpi must be at least 72")
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    factories = (
        ("family_actionability_heatmap", create_family_actionability_heatmap),
        ("gene_family_landscape", create_gene_family_landscape),
        ("family_actionability_stack", create_family_actionability_stack),
        ("rare_disease_family_heatmap", create_rare_disease_family_heatmap),
    )
    plt, _patch = _matplotlib()
    paths: list[Path] = []
    with tempfile.TemporaryDirectory(prefix=".global-plots-", dir=destination) as temporary:
        temp_dir = Path(temporary)
        for stem, factory in factories:
            figure = factory(report)
            try:
                for file_format in normalized_formats:
                    temporary_path = temp_dir / f"{stem}.{file_format}"
                    figure.savefig(temporary_path, dpi=dpi, bbox_inches="tight")
                    final_path = destination / temporary_path.name
                    temporary_path.replace(final_path)
                    paths.append(final_path)
            finally:
                plt.close(figure)
    return tuple(paths)


def _empty_panel(ax: Any, title: str) -> None:
    """Role: Render an empty-state panel placeholder."""
    ax.set_title(title, loc="left", fontweight="bold")
    ax.text(0.5, 0.5, "No compatible cached analyses", ha="center", va="center", transform=ax.transAxes)
    ax.set_axis_off()
