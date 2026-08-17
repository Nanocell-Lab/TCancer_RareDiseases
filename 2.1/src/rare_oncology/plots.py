"""Optional Matplotlib visualizations for a gene analysis."""

from __future__ import annotations

import sys
import tempfile
import warnings
import math
from functools import lru_cache
from pathlib import Path
from textwrap import fill
from typing import Any

from .models import AnalysisResult
from .plot_data import (
    heatmap_data,
    LollipopResidue,
    lollipop_label_residues,
    lollipop_residues,
    RareDiseaseHeatmapData,
    rare_disease_heatmap_data,
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
from .string_network import StringEdge, StringNetwork


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
            "Plots require Matplotlib. Install it with: "
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
            "The lollipop plot requires Geneview. Install it with: "
            "python -m pip install -e '.[plots]'"
        ) from exc
    return pd, lolliplot, GenomicInterval


def create_lollipop_plot(
    result: AnalysisResult,
    *,
    max_mutations: int = 80,
    max_labels: int = 15,
    include_synonymous: bool = False,
    protein_annotation: ProteinAnnotation | None = None,
) -> Any:
    """Role: Create the lollipop plot visualization or result object."""
    _validate_maximum(max_mutations, "max_mutations")
    _validate_nonnegative(max_labels, "max_labels")
    annotation = _resolve_protein_annotation(result, protein_annotation)
    plt, _line2d, _circle = _matplotlib()
    fig, ax = plt.subplots(figsize=(13, 6), layout="constrained")
    try:
        _draw_lollipop(
            ax,
            result,
            annotation,
            max_mutations=max_mutations,
            max_labels=max_labels,
            include_synonymous=include_synonymous,
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


def create_rare_disease_heatmap_plot(
    result: AnalysisResult, *, max_mutations: int = 30
) -> Any:
    """Role: Create the rare disease heatmap plot visualization or result object."""
    _validate_rare_disease_maximum(max_mutations)
    data = rare_disease_heatmap_data(result, max_mutations=max_mutations)
    plt, _line2d, _circle = _matplotlib()
    width = max(12.0, min(20.0, 5.0 + len(data.mutations) * 0.45))
    height = max(6.0, min(24.0, 3.5 + len(data.diseases) * 0.55))
    fig, ax = plt.subplots(figsize=(width, height), layout="constrained")
    try:
        image = _draw_rare_disease_heatmap(
            ax,
            result,
            max_mutations=max_mutations,
            data=data,
        )
        if image is not None:
            fig.colorbar(
                image,
                ax=ax,
                label="ClinVar association count",
            )
    except Exception:
        plt.close(fig)
        raise
    return fig


def create_string_network_plot(network: StringNetwork) -> Any:
    """Create a deterministic, dependency-light STRING association graph."""

    plt, _line2d, _circle = _matplotlib()
    fig, ax = plt.subplots(figsize=(12, 10), layout="constrained")
    try:
        _draw_string_network(ax, network)
    except Exception:
        plt.close(fig)
        raise
    return fig


def generate_string_network_plots(
    network: StringNetwork,
    output_dir: str | Path,
    *,
    formats: tuple[str, ...] = ("png", "svg"),
    dpi: int = 200,
) -> tuple[Path, ...]:
    """Role: Generate and persist string network plots outputs."""
    normalized_formats = tuple(dict.fromkeys(value.lower() for value in formats))
    unsupported = set(normalized_formats) - {"png", "svg"}
    if not normalized_formats or unsupported:
        raise ValueError("formats must contain only 'png' and/or 'svg'")
    if dpi < 72:
        raise ValueError("dpi must be at least 72")
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    figure = create_string_network_plot(network)
    plt, _line2d, _circle = _matplotlib()
    written: list[Path] = []
    try:
        with tempfile.TemporaryDirectory(
            prefix=".rare-oncology-string-", dir=destination
        ) as tmp:
            temporary = Path(tmp)
            staged: list[tuple[Path, Path]] = []
            for image_format in normalized_formats:
                filename = f"{network.gene.lower()}_string_network.{image_format}"
                temporary_path = temporary / filename
                final_path = destination / filename
                figure.savefig(
                    temporary_path,
                    dpi=dpi,
                    bbox_inches="tight",
                    metadata={
                        "Creator": "rare-oncology 2.1",
                        "Title": f"{network.gene} STRING network v{network.string_version}",
                    },
                )
                staged.append((temporary_path, final_path))
            for temporary_path, final_path in staged:
                temporary_path.replace(final_path)
                written.append(final_path)
    finally:
        plt.close(figure)
    return tuple(written)


def _draw_string_network(ax: Any, network: StringNetwork) -> None:
    """Role: Draw string network on the supplied plotting axis."""
    _plt, Line2D, _circle = _matplotlib()
    nodes = {node.string_id: node for node in network.nodes}
    query_id = network.mapped_string_id
    if query_id not in nodes:
        raise ValueError("The STRING network does not contain the query node")
    adjacency: dict[str, set[str]] = {node_id: set() for node_id in nodes}
    degree: dict[str, int] = {node_id: 0 for node_id in nodes}
    for edge in network.edges:
        adjacency.setdefault(edge.string_id_a, set()).add(edge.string_id_b)
        adjacency.setdefault(edge.string_id_b, set()).add(edge.string_id_a)
        degree[edge.string_id_a] = degree.get(edge.string_id_a, 0) + 1
        degree[edge.string_id_b] = degree.get(edge.string_id_b, 0) + 1

    direct = sorted(
        adjacency.get(query_id, set()),
        key=lambda node_id: nodes[node_id].preferred_name,
    )
    secondary = sorted(
        (node_id for node_id in nodes if node_id != query_id and node_id not in direct),
        key=lambda node_id: nodes[node_id].preferred_name,
    )
    positions: dict[str, tuple[float, float]] = {query_id: (0.0, 0.0)}
    positions.update(_ring_positions(direct, radius=1.0, offset=math.pi / 2))
    positions.update(_ring_positions(secondary, radius=1.65, offset=0.0))

    displayed_edges = _string_display_edges(network)
    for edge in sorted(displayed_edges, key=lambda item: item.score):
        start = positions.get(edge.string_id_a)
        end = positions.get(edge.string_id_b)
        if start is None or end is None:
            continue
        score = edge.score
        query_edge = query_id in {edge.string_id_a, edge.string_id_b}
        ax.plot(
            (start[0], end[0]),
            (start[1], end[1]),
            color=(
                (0.92, 0.35, 0.05, 0.28 + score * 0.64)
                if query_edge
                else (0.39, 0.52, 0.65, 0.12 + score * 0.32)
            ),
            linewidth=(0.8 + score * 2.4 if query_edge else 0.4 + score),
            zorder=2 if query_edge else 1,
        )

    partner_ids = [node_id for node_id in positions if node_id != query_id]
    if partner_ids:
        ax.scatter(
            [positions[node_id][0] for node_id in partner_ids],
            [positions[node_id][1] for node_id in partner_ids],
            s=[140 + degree.get(node_id, 0) * 28 for node_id in partner_ids],
            c="#60A5FA",
            edgecolors="#1E3A8A",
            linewidths=1.0,
            zorder=3,
        )
    ax.scatter(
        [0.0],
        [0.0],
        s=[320 + degree.get(query_id, 0) * 32],
        c="#F97316",
        edgecolors="#9A3412",
        linewidths=1.8,
        zorder=4,
    )
    for node_id, (x, y) in positions.items():
        node = nodes[node_id]
        ax.text(
            x,
            y + (0.08 if node_id == query_id else 0.065),
            node.preferred_name,
            ha="center",
            va="bottom",
            fontsize=10 if node_id == query_id else 8,
            fontweight="bold" if node_id == query_id else "normal",
            color="#7C2D12" if node_id == query_id else "#0F172A",
            bbox={
                "boxstyle": "round,pad=0.18",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.78,
            },
            zorder=5,
        )

    network_label = (
        "functional" if network.network_type == "functional" else "physical"
    )
    ax.set_title(
        f"{network.gene}: {network_label} protein associations "
        f"in STRING v{network.string_version}\n"
        f"Homo sapiens (9606) · minimum score {network.required_score}/1000 · "
        f"{len(network.nodes)} nodes · {len(network.edges)} edges",
        fontsize=14,
        fontweight="bold",
        pad=18,
    )
    ax.text(
        0.01,
        0.01,
        f"Showing {len(displayed_edges)} of {len(network.edges)} edges: "
        "all query-gene edges plus the strongest secondary edges. "
        "Association ≠ physical binding or causality.",
        transform=ax.transAxes,
        fontsize=8,
        color="#475569",
        ha="left",
        va="bottom",
    )
    ax.legend(
        handles=[
            Line2D(
                [0], [0], marker="o", color="none", markerfacecolor="#F97316",
                markeredgecolor="#9A3412", markersize=10, label="Query gene"
            ),
            Line2D(
                [0], [0], marker="o", color="none", markerfacecolor="#60A5FA",
                markeredgecolor="#1E3A8A", markersize=9, label="Associated protein"
            ),
            Line2D(
                [0], [0], color="#EA580C", linewidth=3,
                label="Association with the query gene"
            ),
            Line2D(
                [0], [0], color="#64748B", linewidth=1.5,
                label="Strong secondary association"
            ),
        ],
        loc="upper right",
        fontsize=8,
    )
    limit = 1.95 if secondary else 1.28
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal")
    ax.axis("off")


def _ring_positions(
    node_ids: list[str], *, radius: float, offset: float
) -> dict[str, tuple[float, float]]:
    """Role: Compute ring positions for the plots workflow."""
    if not node_ids:
        return {}
    return {
        node_id: (
            radius * math.cos(offset + 2 * math.pi * index / len(node_ids)),
            radius * math.sin(offset + 2 * math.pi * index / len(node_ids)),
        )
        for index, node_id in enumerate(node_ids)
    }


def _string_display_edges(network: StringNetwork) -> tuple[StringEdge, ...]:
    """Role: Compute string display edges for the plots workflow."""
    query_id = network.mapped_string_id
    query_edges = [
        edge
        for edge in network.edges
        if query_id in {edge.string_id_a, edge.string_id_b}
    ]
    secondary = [
        edge
        for edge in network.edges
        if query_id not in {edge.string_id_a, edge.string_id_b}
    ]
    secondary_limit = max(8, (len(network.nodes) - 1) // 2)
    secondary.sort(
        key=lambda edge: (
            -edge.score,
            edge.preferred_name_a,
            edge.preferred_name_b,
        )
    )
    return tuple(query_edges + secondary[:secondary_limit])


def create_dashboard(
    result: AnalysisResult,
    *,
    max_lollipop_mutations: int = 80,
    max_lollipop_labels: int = 15,
    max_heatmap_mutations: int = 30,
    include_synonymous_lollipop: bool = False,
    protein_annotation: ProteinAnnotation | None = None,
) -> Any:
    """Role: Create the dashboard visualization or result object."""
    _validate_maximum(max_lollipop_mutations, "max_lollipop_mutations")
    _validate_nonnegative(max_lollipop_labels, "max_lollipop_labels")
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
            max_labels=min(max_lollipop_labels, 8),
            include_synonymous=include_synonymous_lollipop,
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
    max_lollipop_mutations: int = 80,
    max_lollipop_labels: int = 15,
    max_heatmap_mutations: int = 30,
    max_rare_disease_mutations: int = 30,
    include_synonymous_lollipop: bool = False,
    dpi: int = 200,
    protein_annotation: ProteinAnnotation | None = None,
) -> dict[str, tuple[Path, ...]]:
    """Create five plots plus a dashboard and return all written paths."""

    normalized_formats = tuple(dict.fromkeys(value.lower() for value in formats))
    unsupported = set(normalized_formats) - {"png", "svg"}
    if not normalized_formats or unsupported:
        raise ValueError("formats must contain only 'png' and/or 'svg'")
    if dpi < 72:
        raise ValueError("dpi must be at least 72")
    _validate_maximum(max_lollipop_mutations, "max_lollipop_mutations")
    _validate_nonnegative(max_lollipop_labels, "max_lollipop_labels")
    _validate_maximum(max_heatmap_mutations, "max_heatmap_mutations")
    _validate_rare_disease_maximum(max_rare_disease_mutations)
    annotation = _resolve_protein_annotation(result, protein_annotation)

    creators = {
        "lollipop": lambda: create_lollipop_plot(
            result,
            max_mutations=max_lollipop_mutations,
            max_labels=max_lollipop_labels,
            include_synonymous=include_synonymous_lollipop,
            protein_annotation=annotation,
        ),
        "venn": lambda: create_venn_plot(result),
        "recurrence": lambda: create_recurrence_plot(result),
        "heatmap": lambda: create_heatmap_plot(
            result, max_mutations=max_heatmap_mutations
        ),
        "rare_disease_heatmap": lambda: create_rare_disease_heatmap_plot(
            result, max_mutations=max_rare_disease_mutations
        ),
        "dashboard": lambda: create_dashboard(
            result,
            max_lollipop_mutations=max_lollipop_mutations,
            max_lollipop_labels=max_lollipop_labels,
            max_heatmap_mutations=max_heatmap_mutations,
            include_synonymous_lollipop=include_synonymous_lollipop,
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
                            "Creator": "rare-oncology 2.1",
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
    max_labels: int,
    include_synonymous: bool,
) -> None:
    """Role: Draw lollipop on the supplied plotting axis."""
    canonical_limit = annotation.length if annotation.accession else None
    residues = lollipop_residues(
        result,
        max_residues=max_mutations,
        max_position=canonical_limit,
        include_synonymous=include_synonymous,
    )
    eligible_records = tuple(
        record
        for record in result.mutation_plot_data
        if record.amino_acid_position is not None
        and (canonical_limit is None or record.amino_acid_position <= canonical_limit)
        and (include_synonymous or not record.protein_change.endswith("="))
    )
    total_residues = len(
        {record.amino_acid_position for record in eligible_records}
    )
    out_of_range = sum(
        record.amino_acid_position is not None
        and canonical_limit is not None
        and record.amino_acid_position > canonical_limit
        for record in result.mutation_plot_data
    )
    synonymous_omitted = sum(
        record.amino_acid_position is not None
        and (canonical_limit is None or record.amino_acid_position <= canonical_limit)
        and record.protein_change.endswith("=")
        for record in result.mutation_plot_data
    ) if not include_synonymous else 0
    suffix = (
        f" ({len(residues)} residues from {total_residues}; "
        f"{sum(len(residue.protein_changes) for residue in residues)} changes)"
        if total_residues > len(residues)
        else ""
    )
    provenance = (
        f" · UniProt {annotation.accession}"
        if annotation.accession
        else " · inferred length (no UniProt)"
    )
    range_note = (
        f" · {out_of_range} out-of-range change"
        f"{'s' if out_of_range != 1 else ''} omitted"
        if out_of_range
        else ""
    )
    synonymous_note = (
        f" · {synonymous_omitted} synonymous changes omitted"
        if synonymous_omitted
        else ""
    )
    title = f"{result.gene}: variant protein map{suffix}{provenance}"
    notes = [note.removeprefix(" · ") for note in (range_note, synonymous_note) if note]
    if notes:
        title += "\n" + " · ".join(notes)
    if not residues:
        ax.set_title(title)
        _empty_axis(ax, "No positionable protein changes")
        return
    pd, lolliplot, GenomicInterval = _geneview()
    maximum_position = max(residue.position for residue in residues)
    protein_length = (
        annotation.length
        if annotation.accession
        else max(annotation.length, maximum_position)
    )
    variants = pd.DataFrame(
        [
            {
                "chrom": result.gene,
                "start": residue.position,
                "end": residue.position + 1,
                "score": max(residue.cancer_patient_count, 1),
                "label": "",
                "pie_values": [
                    residue.rare_count,
                    residue.oncology_count,
                    residue.shared_count,
                ],
                "pie_colors": [
                    COLORS["rare"],
                    COLORS["oncology"],
                    COLORS["shared"],
                ],
            }
            for residue in residues
        ]
    )
    domains = _geneview_domains(pd, result.gene, annotation.domains)
    lolliplot(
        variants,
        features=domains,
        type="pie.stack",
        region=GenomicInterval(result.gene, 1, protein_length + 1),
        ax=ax,
        show_yaxis=False,
        label_on_feature=True,
        node_size=0.012,
        cex=1.0,
        jitter=None,
    )
    _soften_geneview_stems(ax)
    _draw_exact_residue_labels(
        ax,
        lollipop_label_residues(residues, max_labels=max_labels),
        protein_length,
    )
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
    ax.set_xlabel("Protein change")
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
        f" (top {len(data.mutations)} of {total_oncology})"
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
    ax.set_ylabel("Protein change")
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


def _draw_rare_disease_heatmap(
    ax: Any,
    result: AnalysisResult,
    *,
    max_mutations: int,
    data: RareDiseaseHeatmapData | None = None,
) -> Any | None:
    """Role: Draw rare disease heatmap on the supplied plotting axis."""
    if data is None:
        data = rare_disease_heatmap_data(result, max_mutations=max_mutations)
    total_rare = sum(
        record.in_rare_disease for record in result.mutation_plot_data
    )
    suffix = (
        f" (top {len(data.mutations)} of {total_rare})"
        if total_rare > len(data.mutations)
        else ""
    )
    ax.set_title(
        f"{result.gene}: mutations by rare disease{suffix}\n"
        "[C] indicates a mutation shared with cancer",
    )
    if not data.mutations or not data.diseases:
        _empty_axis(
            ax,
            "No rare-disease associations are available for the heatmap",
        )
        return None
    image = ax.imshow(
        data.values,
        cmap="YlGnBu",
        interpolation="nearest",
        aspect="auto",
        vmin=0,
    )
    mutation_labels = [
        f"[C] p.{mutation}" if shared else f"p.{mutation}"
        for mutation, shared in zip(
            data.mutations, data.shared_mutations, strict=True
        )
    ]
    ax.set_xticks(
        range(len(data.mutations)),
        labels=mutation_labels,
        rotation=55,
        ha="right",
    )
    for tick, shared in zip(
        ax.get_xticklabels(), data.shared_mutations, strict=True
    ):
        if shared:
            tick.set_color(COLORS["shared"])
            tick.set_fontweight("bold")
    ax.set_yticks(
        range(len(data.diseases)),
        labels=[
            fill(disease, width=38, break_long_words=False)
            for disease in data.diseases
        ],
    )
    ax.set_xlabel("Protein change; [C] indicates a cancer match")
    ax.set_ylabel("Associated rare disease in ClinVar")
    if len(data.mutations) * len(data.diseases) <= 240:
        maximum = max(max(row) for row in data.values) or 1
        for row_index, row in enumerate(data.values):
            for column_index, value in enumerate(row):
                if value == 0:
                    continue
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
        raise ValueError(f"{name} must be positive")


def _validate_nonnegative(value: int, name: str) -> None:
    """Role: Validate nonnegative constraints before execution."""
    if value < 0:
        raise ValueError(f"{name} cannot be negative")


def _validate_rare_disease_maximum(value: int) -> None:
    """Role: Validate rare disease maximum constraints before execution."""
    if not 1 <= value <= 30:
        raise ValueError("max_rare_disease_mutations must be between 1 and 30")


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
            f"UniProt could not be queried for {result.gene}: {exc}",
            RuntimeWarning,
            stacklevel=2,
        )
        annotation = None
    if annotation is None and not lookup_failed:
        warnings.warn(
            f"No reviewed human UniProt entry is available for {result.gene}; "
            "an inferred length without domains or secondary structure will be used.",
            RuntimeWarning,
            stacklevel=2,
        )
    return annotation or inferred_annotation(result.gene, positions)


def _geneview_domains(pd: Any, gene: str, domains: tuple[ProteinFeature, ...]) -> Any:
    """Role: Load the optional Geneview plotting dependencies."""
    if not domains:
        return None
    palette = {
        "signal": "#FDE68A",
        "topology": "#DBEAFE",
        "transmembrane": "#FDBA74",
        "domain": "#C4B5FD",
    }
    return pd.DataFrame(
        [
            {
                "chrom": gene,
                "start": domain.start,
                "end": domain.end + 1,
                "name": domain.label,
                "fill": palette.get(domain.kind, "#A7F3D0"),
                "color": "#475569",
                "height": 0.065,
                "feature": "domain",
                "feature_layer_id": (
                    "topology" if domain.kind == "topology" else "functional"
                ),
            }
            for domain in domains
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
    for y, label in (
        (0.48, "Variants"),
        (0.08, "Domains"),
        (-0.16, "Sec. structure"),
    ):
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


def _draw_exact_residue_labels(
    ax: Any,
    residues: tuple[LollipopResidue, ...],
    protein_length: int,
) -> None:
    """Place a limited number of labels without changing residue coordinates."""

    if not residues:
        return
    levels = (0.86, 0.94, 1.02, 1.10)
    last_right = [float("-inf")] * len(levels)
    axes_width_px = max(ax.get_window_extent().width, 1.0)
    data_per_pixel = protein_length / axes_width_px
    for residue in sorted(residues, key=lambda item: item.position):
        x = residue.position
        extra = len(residue.protein_changes) - 1
        text = f"p.{residue.label_change}"
        if extra:
            text += f" (+{extra})"
        if residue.cancer_patient_count:
            text += f" · n={residue.cancer_patient_count}"
        half_width = max(18.0, len(text) * 3.2 * data_per_pixel)
        level_index = next(
            (
                index
                for index, previous_right in enumerate(last_right)
                if x - half_width >= previous_right + 8.0
            ),
            min(range(len(levels)), key=last_right.__getitem__),
        )
        label_y = levels[level_index]
        ax.vlines(
            x,
            0.72,
            label_y - 0.015,
            color="#94A3B8",
            linewidth=0.7,
            linestyle=":",
            zorder=4,
        )
        ax.text(
            x,
            label_y,
            text,
            ha="center",
            va="bottom",
            fontsize=7,
            color="#111827",
            clip_on=False,
            zorder=7,
        )
        last_right[level_index] = x + half_width


def _soften_geneview_stems(ax: Any) -> None:
    """Reduce the visual weight of Geneview's exact-position guide stems."""

    for line in ax.lines:
        xdata = line.get_xdata()
        if len(xdata) != 2 or float(xdata[0]) != float(xdata[1]):
            continue
        line.set_color("#CBD5E1")
        line.set_linewidth(0.65)
        line.set_alpha(0.85)


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
