"""Plotly figures for the interactive dashboard."""

from __future__ import annotations

import math
from collections import Counter
from typing import Any, Mapping, Sequence

import plotly.graph_objects as go


COLORS = {
    "Rare disease only": "#3B82F6",
    "Shared": "#F97316",
    "Oncology only": "#A855F7",
}


def empty_figure(message: str) -> go.Figure:
    """Role: Provide the empty figure operation for the figures workflow."""
    figure = go.Figure()
    figure.add_annotation(
        text=message,
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font={"size": 16, "color": "#64748B"},
    )
    figure.update_xaxes(visible=False)
    figure.update_yaxes(visible=False)
    figure.update_layout(height=320, margin={"l": 20, "r": 20, "t": 40, "b": 20})
    return figure


def overlap_figure(rows: Sequence[Mapping[str, Any]]) -> go.Figure:
    """Role: Provide the overlap figure operation for the figures workflow."""
    figure = go.Figure(
        go.Bar(
            x=[row["category"] for row in rows],
            y=[row["count"] for row in rows],
            marker_color=[COLORS[str(row["category"])] for row in rows],
            text=[row["count"] for row in rows],
            textposition="outside",
            hovertemplate="%{x}: %{y}<extra></extra>",
        )
    )
    figure.update_layout(
        title="Mutation-set overlap",
        xaxis_title=None,
        yaxis_title="Unique protein changes",
        showlegend=False,
        height=360,
        margin={"l": 50, "r": 20, "t": 60, "b": 70},
    )
    return figure


def lollipop_figure(rows: Sequence[Mapping[str, Any]], gene: str) -> go.Figure:
    """Role: Provide the lollipop figure operation for the figures workflow."""
    if not rows:
        return empty_figure("No positioned protein changes are available.")
    figure = go.Figure()
    max_position = max(int(row["position"]) for row in rows)
    figure.add_trace(
        go.Scatter(
            x=[1, max_position],
            y=[0, 0],
            mode="lines",
            line={"color": "#64748B", "width": 7},
            hoverinfo="skip",
            showlegend=False,
        )
    )
    for category in ("Shared", "Oncology only", "Rare disease only"):
        selected = [row for row in rows if row["category"] == category]
        if not selected:
            continue
        line_x: list[int | None] = []
        line_y: list[float | None] = []
        for row in selected:
            line_x.extend([int(row["position"]), int(row["position"]), None])
            line_y.extend([0.0, float(row["height"]), None])
        figure.add_trace(
            go.Scatter(
                x=line_x,
                y=line_y,
                mode="lines",
                line={"color": COLORS[category], "width": 2},
                hoverinfo="skip",
                showlegend=False,
            )
        )
        figure.add_trace(
            go.Scatter(
                x=[row["position"] for row in selected],
                y=[row["height"] for row in selected],
                mode="markers",
                name=category,
                marker={
                    "color": COLORS[category],
                    "size": [
                        min(24, 10 + 2 * int(row["mutation_count"])) for row in selected
                    ],
                    "line": {
                        "color": ["#111827" if row["confirmed"] else "#FFFFFF" for row in selected],
                        "width": [2.5 if row["confirmed"] else 1 for row in selected],
                    },
                },
                customdata=[
                    [
                        row["protein_changes"],
                        row["cancer_patient_count"],
                        row["diseases"] or "—",
                        row["projects"] or "—",
                        "Yes" if row["confirmed"] else "No",
                    ]
                    for row in selected
                ],
                hovertemplate=(
                    "Residue %{x}<br>Variants: %{customdata[0]}"
                    "<br>Cancer patients: %{customdata[1]}"
                    "<br>Rare diseases: %{customdata[2]}"
                    "<br>Projects: %{customdata[3]}"
                    "<br>Genomic confirmation: %{customdata[4]}<extra>%{fullData.name}</extra>"
                ),
            )
        )
    figure.update_layout(
        title=f"{gene} interactive protein lollipop",
        xaxis_title="Residue position",
        yaxis_title="Recurrence (log-scaled stem height)",
        hovermode="closest",
        height=520,
        legend={"orientation": "h", "y": 1.12, "x": 0},
        margin={"l": 65, "r": 25, "t": 90, "b": 60},
    )
    figure.update_yaxes(showticklabels=False, zeroline=False)
    return figure


def recurrence_figure(rows: Sequence[Mapping[str, Any]], *, top: int = 30) -> go.Figure:
    """Role: Provide the recurrence figure operation for the figures workflow."""
    selected = list(rows[:top])
    if not selected:
        return empty_figure("No shared recurrent mutations are available.")
    selected.reverse()
    figure = go.Figure(
        go.Bar(
            x=[row["cancer_patients"] for row in selected],
            y=[f"p.{row['protein_change']}" for row in selected],
            orientation="h",
            marker_color=[
                "#EA580C"
                if row["match_level"] == "confirmed_genomic_and_protein"
                else "#FDBA74"
                for row in selected
            ],
            customdata=[
                [row["rare_diseases"], row["oncology_projects"], row["match_level"]]
                for row in selected
            ],
            hovertemplate=(
                "%{y}<br>Cancer patients: %{x}<br>Rare diseases: %{customdata[0]}"
                "<br>Projects: %{customdata[1]}<br>Match: %{customdata[2]}<extra></extra>"
            ),
        )
    )
    figure.update_layout(
        title=f"Top {len(selected)} shared mutations by cancer patient count",
        xaxis_title="Unique cancer patients",
        yaxis_title=None,
        height=max(420, 24 * len(selected)),
        margin={"l": 90, "r": 25, "t": 60, "b": 50},
    )
    return figure


def heatmap_figure(
    matrix: Mapping[str, Any],
    *,
    row_key: str,
    column_key: str,
    title: str,
    colorbar_title: str,
) -> go.Figure:
    """Role: Provide the heatmap figure operation for the figures workflow."""
    rows = list(matrix.get(row_key, []))
    columns = list(matrix.get(column_key, []))
    values = list(matrix.get("values", []))
    if not rows or not columns:
        return empty_figure("Not enough data to build this heatmap.")
    shared = list(matrix.get("shared", []))
    if row_key == "mutations":
        y_labels = [
            f"{'★ ' if index < len(shared) and shared[index] else ''}p.{label}"
            for index, label in enumerate(rows)
        ]
        x_labels = columns
        z = values
    else:
        y_labels = rows
        x_labels = [
            f"{'★ ' if index < len(shared) and shared[index] else ''}p.{label}"
            for index, label in enumerate(columns)
        ]
        z = values
    figure = go.Figure(
        go.Heatmap(
            z=z,
            x=x_labels,
            y=y_labels,
            colorscale="Blues",
            xgap=1,
            ygap=1,
            colorbar={"title": colorbar_title},
            hovertemplate="%{y}<br>%{x}<br>Count: %{z}<extra></extra>",
        )
    )
    figure.update_layout(
        title=title,
        height=max(480, min(1000, 25 * len(y_labels) + 180)),
        margin={"l": 160 if row_key == "diseases" else 90, "r": 30, "t": 70, "b": 130},
    )
    figure.update_xaxes(tickangle=-45)
    return figure


def alphamissense_figure(rows: Sequence[Mapping[str, Any]], gene: str) -> go.Figure:
    """Role: Provide the alphamissense figure operation for the figures workflow."""
    positioned = [row for row in rows if row.get("position") is not None and row.get("score") is not None]
    if not positioned:
        return empty_figure("No AlphaMissense scores are stored for this snapshot.")
    predictions = sorted({str(row["prediction"]) for row in positioned})
    palette = ["#DC2626", "#F59E0B", "#2563EB", "#64748B", "#16A34A"]
    figure = go.Figure()
    for prediction, color in zip(predictions, palette, strict=False):
        selected = [row for row in positioned if row["prediction"] == prediction]
        figure.add_trace(
            go.Scatter(
                x=[row["position"] for row in selected],
                y=[row["score"] for row in selected],
                mode="markers",
                name=prediction,
                marker={
                    "color": color,
                    "size": [13 if row["shared_with_oncology"] else 8 for row in selected],
                    "symbol": ["diamond" if row["shared_with_oncology"] else "circle" for row in selected],
                    "opacity": 0.82,
                },
                text=[f"p.{row['protein_change']}" for row in selected],
                customdata=[["Yes" if row["shared_with_oncology"] else "No"] for row in selected],
                hovertemplate=(
                    "%{text}<br>Score: %{y:.3f}<br>Shared with oncology: %{customdata[0]}"
                    "<extra>%{fullData.name}</extra>"
                ),
            )
        )
    figure.add_hrect(y0=0.56, y1=1.0, fillcolor="#FEE2E2", opacity=0.18, line_width=0)
    figure.update_layout(
        title=f"{gene} AlphaMissense predictions",
        xaxis_title="Residue position",
        yaxis_title="AlphaMissense pathogenicity score",
        yaxis_range=[0, 1],
        height=480,
        legend={"orientation": "h", "y": 1.15},
        margin={"l": 65, "r": 25, "t": 90, "b": 60},
    )
    return figure


def therapy_summary_figure(rows: Sequence[Mapping[str, Any]]) -> go.Figure:
    """Role: Provide the therapy summary figure operation for the figures workflow."""
    if not rows:
        return empty_figure("No therapy snapshot is stored for this gene.")
    counts = Counter((str(row["layer"]), str(row.get("status") or "Unspecified")) for row in rows)
    layers = list(dict.fromkeys(layer for layer, _status in counts))
    statuses = sorted({status for _layer, status in counts})
    figure = go.Figure()
    for status in statuses:
        figure.add_trace(
            go.Bar(
                name=status,
                x=layers,
                y=[counts.get((layer, status), 0) for layer in layers],
                hovertemplate="%{x}<br>%{fullData.name}: %{y}<extra></extra>",
            )
        )
    figure.update_layout(
        title="Therapy associations by evidence layer and status",
        barmode="stack",
        yaxis_title="Associations",
        xaxis_title=None,
        height=440,
        margin={"l": 60, "r": 25, "t": 60, "b": 100},
        legend={"orientation": "h", "y": 1.18},
    )
    return figure


def string_network_figure(graph: Mapping[str, Any], gene: str) -> go.Figure:
    """Role: Provide the string network figure operation for the figures workflow."""
    nodes = list(graph.get("nodes", []))
    edges = list(graph.get("edges", []))
    query_id = str(graph.get("query_id") or "")
    if not nodes or not query_id:
        return empty_figure("No stored STRING network is available for this gene.")

    query = next(
        (node for node in nodes if str(node.get("string_id")) == query_id),
        nodes[0],
    )
    partners = sorted(
        (node for node in nodes if node is not query),
        key=lambda node: str(node.get("preferred_name") or node.get("string_id")),
    )
    positions: dict[str, tuple[float, float]] = {query_id: (0.0, 0.0)}
    count = max(1, len(partners))
    for index, node in enumerate(partners):
        angle = (2.0 * math.pi * index / count) + math.pi / 2.0
        radius = 1.0 + 0.09 * (index % 2)
        positions[str(node.get("string_id"))] = (
            radius * math.cos(angle),
            radius * math.sin(angle),
        )

    figure = go.Figure()
    for direct, color, width, name in (
        (False, "rgba(100,116,139,.34)", 1.2, "Partner–partner"),
        (True, "rgba(249,115,22,.70)", 2.3, f"Direct to {gene}"),
    ):
        x_values: list[float | None] = []
        y_values: list[float | None] = []
        for edge in edges:
            node_a = str(edge.get("string_id_a") or "")
            node_b = str(edge.get("string_id_b") or "")
            is_direct = query_id in {node_a, node_b}
            if is_direct != direct or node_a not in positions or node_b not in positions:
                continue
            x_values.extend([positions[node_a][0], positions[node_b][0], None])
            y_values.extend([positions[node_a][1], positions[node_b][1], None])
        if x_values:
            figure_trace = go.Scatter(
                x=x_values,
                y=y_values,
                mode="lines",
                name=name,
                line={"color": color, "width": width},
                hoverinfo="skip",
            )
            figure.add_trace(figure_trace)

    for is_query, selected, color, size, label in (
        (False, partners, "#2563EB", 15, "STRING partner"),
        (True, [query], "#F97316", 27, gene),
    ):
        if not selected:
            continue
        figure.add_trace(
            go.Scatter(
                x=[positions[str(node.get("string_id"))][0] for node in selected],
                y=[positions[str(node.get("string_id"))][1] for node in selected],
                mode="markers+text",
                name=label,
                text=[str(node.get("preferred_name") or "") for node in selected],
                textposition="top center",
                marker={
                    "color": color,
                    "size": [
                        size if is_query else min(24, size + int(node.get("degree", 0)))
                        for node in selected
                    ],
                    "line": {"color": "white", "width": 1.5},
                },
                customdata=[
                    [node.get("string_id"), node.get("degree", 0)] for node in selected
                ],
                hovertemplate=(
                    "%{text}<br>STRING ID: %{customdata[0]}"
                    "<br>Network degree: %{customdata[1]}<extra>%{fullData.name}</extra>"
                ),
            )
        )
    figure.update_layout(
        title=f"{gene} functional partners in STRING",
        height=650,
        hovermode="closest",
        xaxis={"visible": False},
        yaxis={"visible": False, "scaleanchor": "x", "scaleratio": 1},
        legend={"orientation": "h", "y": 1.08, "x": 0},
        margin={"l": 30, "r": 30, "t": 85, "b": 35},
        plot_bgcolor="rgba(255,255,255,.7)",
    )
    return figure


def string_partner_scores_figure(
    rows: Sequence[Mapping[str, Any]], *, top: int = 20
) -> go.Figure:
    """Role: Provide the string partner scores figure operation for the figures workflow."""
    selected = list(rows[:top])
    if not selected:
        return empty_figure("No direct STRING partners are stored.")
    selected.reverse()
    figure = go.Figure(
        go.Bar(
            x=[row["combined_score"] for row in selected],
            y=[row["partner"] for row in selected],
            orientation="h",
            marker_color="#2563EB",
            customdata=[
                [
                    row["experimental_score"],
                    row["database_score"],
                    row["coexpression_score"],
                    row["textmining_score"],
                ]
                for row in selected
            ],
            hovertemplate=(
                "%{y}<br>Combined: %{x:.3f}<br>Experimental: %{customdata[0]:.3f}"
                "<br>Database: %{customdata[1]:.3f}"
                "<br>Coexpression: %{customdata[2]:.3f}"
                "<br>Text mining: %{customdata[3]:.3f}<extra></extra>"
            ),
        )
    )
    figure.update_layout(
        title="Direct partners ranked by combined score",
        xaxis_title="STRING combined score",
        xaxis_range=[0, 1],
        yaxis_title=None,
        height=max(420, 26 * len(selected) + 150),
        margin={"l": 100, "r": 25, "t": 65, "b": 55},
    )
    return figure


def string_evidence_heatmap_figure(
    rows: Sequence[Mapping[str, Any]], *, top: int = 20
) -> go.Figure:
    """Role: Provide the string evidence heatmap figure operation for the figures workflow."""
    selected = list(rows[:top])
    if not selected:
        return empty_figure("No STRING evidence channels are stored.")
    channels = (
        ("experimental_score", "Experiments"),
        ("database_score", "Databases"),
        ("coexpression_score", "Coexpression"),
        ("textmining_score", "Text mining"),
        ("neighborhood_score", "Neighborhood"),
        ("fusion_score", "Fusion"),
        ("phylogenetic_score", "Co-occurrence"),
    )
    figure = go.Figure(
        go.Heatmap(
            z=[[row[key] for key, _label in channels] for row in selected],
            x=[label for _key, label in channels],
            y=[row["partner"] for row in selected],
            zmin=0,
            zmax=1,
            colorscale="YlOrRd",
            xgap=1,
            ygap=1,
            colorbar={"title": "Score"},
            hovertemplate="%{y}<br>%{x}: %{z:.3f}<extra></extra>",
        )
    )
    figure.update_layout(
        title="Evidence channels behind each association",
        height=max(440, 25 * len(selected) + 170),
        margin={"l": 100, "r": 25, "t": 65, "b": 85},
    )
    figure.update_xaxes(tickangle=-25)
    return figure


def structure_overlay_figure(comparison: Mapping[str, Any], gene: str) -> go.Figure:
    """Role: Provide the structure overlay figure operation for the figures workflow."""
    wildtype = list(comparison.get("wildtype_coordinates", []))
    variant = list(comparison.get("aligned_variant_coordinates", []))
    residues = list(comparison.get("residues", []))
    if not wildtype or len(wildtype) != len(variant):
        return empty_figure("No aligned C-alpha coordinates are available.")
    labels = [
        f"{row.get('chain', '_')}:{row.get('residue_number', '?')}"
        for row in residues
    ]
    figure = go.Figure()
    for coordinates, name, color in (
        (wildtype, "Wild type", "#2563EB"),
        (variant, "Variant (aligned)", "#F97316"),
    ):
        figure.add_trace(
            go.Scatter3d(
                x=[point[0] for point in coordinates],
                y=[point[1] for point in coordinates],
                z=[point[2] for point in coordinates],
                mode="lines+markers",
                name=name,
                line={"color": color, "width": 5},
                marker={"color": color, "size": 2},
                text=labels,
                hovertemplate="%{text}<extra>%{fullData.name}</extra>",
            )
        )
    figure.update_layout(
        title=f"{gene} Cα structural overlay after Kabsch alignment",
        height=650,
        scene={
            "xaxis_title": "X (Å)",
            "yaxis_title": "Y (Å)",
            "zaxis_title": "Z (Å)",
            "aspectmode": "data",
        },
        legend={"orientation": "h", "y": 1.02},
        margin={"l": 0, "r": 0, "t": 65, "b": 0},
    )
    return figure


def residue_displacement_figure(comparison: Mapping[str, Any]) -> go.Figure:
    """Role: Provide the residue displacement figure operation for the figures workflow."""
    rows = list(comparison.get("residues", []))
    if not rows:
        return empty_figure("No common residues are available.")
    labels = [
        f"{row.get('chain', '_')}:{row.get('residue_number', '?')}"
        f"{row.get('insertion_code', '')}"
        for row in rows
    ]
    colors = ["#DC2626" if row.get("residue_changed") else "#2563EB" for row in rows]
    figure = go.Figure(
        go.Scatter(
            x=labels,
            y=[row.get("ca_displacement_angstrom", 0.0) for row in rows],
            mode="lines+markers",
            line={"color": "#94A3B8", "width": 1.5},
            marker={"color": colors, "size": 6},
            customdata=[
                [row.get("wildtype_residue"), row.get("variant_residue")]
                for row in rows
            ],
            hovertemplate=(
                "%{x}<br>Cα displacement: %{y:.3f} Å"
                "<br>WT: %{customdata[0]} · Variant: %{customdata[1]}<extra></extra>"
            ),
        )
    )
    figure.update_layout(
        title="Per-residue displacement after global alignment",
        xaxis_title="Chain:residue",
        yaxis_title="Cα displacement (Å)",
        height=430,
        margin={"l": 65, "r": 25, "t": 65, "b": 80},
    )
    figure.update_xaxes(nticks=min(30, len(rows)), tickangle=-45)
    return figure
