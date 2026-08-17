"""Pure transformations from 2.1 JSON snapshots to dashboard-ready rows."""

from __future__ import annotations

import math
import re
from collections import defaultdict
from typing import Any, Iterable, Mapping


Payload = Mapping[str, Any]
_POSITION = re.compile(r"(?:p\.)?[^0-9]*(\d+)")
_SIMPLE_MISSENSE = re.compile(r"^[A-Z](\d+)[A-Z]$")


def mutation_category(record: Payload) -> str:
    """Role: Provide the mutation category operation for the transforms workflow."""
    rare = bool(record.get("in_rare_disease"))
    oncology = bool(record.get("in_oncology"))
    if rare and oncology:
        return "Shared"
    if rare:
        return "Rare disease only"
    return "Oncology only"


def overlap_counts(result: Payload) -> list[dict[str, Any]]:
    """Role: Provide the overlap counts operation for the transforms workflow."""
    counts = {"Rare disease only": 0, "Shared": 0, "Oncology only": 0}
    for record in _records(result, "mutation_plot_data"):
        counts[mutation_category(record)] += 1
    return [{"category": category, "count": count} for category, count in counts.items()]


def lollipop_rows(result: Payload, *, max_residues: int = 100) -> list[dict[str, Any]]:
    """Role: Provide the lollipop rows operation for the transforms workflow."""
    if max_residues < 1:
        raise ValueError("max_residues must be positive")
    grouped: dict[int, list[Payload]] = defaultdict(list)
    for record in _records(result, "mutation_plot_data"):
        position = _integer(record.get("amino_acid_position"))
        if position is not None and position > 0:
            grouped[position].append(record)

    rows: list[dict[str, Any]] = []
    for position, records in grouped.items():
        categories = {mutation_category(record) for record in records}
        if "Shared" in categories:
            category = "Shared"
        elif "Oncology only" in categories:
            category = "Oncology only"
        else:
            category = "Rare disease only"
        protein_changes = sorted(
            {str(record.get("protein_change", "")) for record in records if record.get("protein_change")}
        )
        diseases = sorted(
            {
                str(disease)
                for record in records
                for disease in _strings(record.get("diseases"))
            }
        )
        projects = sorted(
            {
                str(project)
                for record in records
                for project, _count in _pairs(record.get("project_patient_counts"))
            }
        )
        patient_count = max(
            (_integer(record.get("cancer_patient_count")) or 0 for record in records),
            default=0,
        )
        height = 1.0 + math.log2(patient_count + 1)
        rows.append(
            {
                "position": position,
                "category": category,
                "protein_changes": ", ".join(protein_changes),
                "mutation_count": len(protein_changes),
                "cancer_patient_count": patient_count,
                "height": height,
                "diseases": "; ".join(diseases),
                "projects": ", ".join(projects),
                "confirmed": any(
                    record.get("match_level") == "confirmed_genomic_and_protein"
                    for record in records
                ),
            }
        )

    priority = {"Shared": 0, "Oncology only": 1, "Rare disease only": 2}
    buckets = {
        category: sorted(
            (row for row in rows if row["category"] == category),
            key=lambda row: (
                not row["confirmed"],
                -int(row["cancer_patient_count"]),
                int(row["position"]),
            ),
        )
        for category in priority
    }
    selected: list[dict[str, Any]] = []
    pattern = ("Shared", "Oncology only", "Shared", "Oncology only", "Rare disease only")
    indices = {category: 0 for category in priority}
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
    return sorted(selected, key=lambda row: int(row["position"]))


def shared_mutation_rows(result: Payload) -> list[dict[str, Any]]:
    """Role: Provide the shared mutation rows operation for the transforms workflow."""
    rows = []
    for record in _records(result, "shared_mutations"):
        rows.append(
            {
                "protein_change": record.get("protein_change"),
                "match_level": record.get("match_level"),
                "cancer_patients": _integer(record.get("cancer_patient_count")) or 0,
                "cancer_events": _integer(record.get("cancer_event_count")) or 0,
                "genomic_matches": _integer(record.get("genomic_match_count")) or 0,
                "rare_diseases": "; ".join(_strings(record.get("diseases"))),
                "orphanet_codes": ", ".join(_strings(record.get("orphanet_codes"))),
                "oncology_projects": ", ".join(_strings(record.get("tcga_projects"))),
                "clinvar_accessions": ", ".join(
                    _strings(record.get("clinvar_accessions"))
                ),
                "classifications": "; ".join(
                    _strings(record.get("germline_classifications"))
                ),
            }
        )
    return sorted(rows, key=lambda row: (-int(row["cancer_patients"]), str(row["protein_change"])))


def oncology_heatmap(result: Payload, *, top: int = 30) -> dict[str, Any]:
    """Role: Provide the oncology heatmap operation for the transforms workflow."""
    if top < 1:
        raise ValueError("top must be positive")
    records = [
        record
        for record in _records(result, "mutation_plot_data")
        if record.get("in_oncology")
    ]
    records.sort(
        key=lambda record: (
            -(_integer(record.get("cancer_patient_count")) or 0),
            str(record.get("protein_change", "")),
        )
    )
    records = records[:top]
    projects = sorted(
        {
            str(project)
            for record in records
            for project, _count in _pairs(record.get("project_patient_counts"))
        }
    )
    mutations = [str(record.get("protein_change", "")) for record in records]
    values = []
    shared = []
    for record in records:
        counts = {str(project): int(count) for project, count in _pairs(record.get("project_patient_counts"))}
        values.append([counts.get(project, 0) for project in projects])
        shared.append(bool(record.get("in_rare_disease")))
    return {"mutations": mutations, "projects": projects, "values": values, "shared": shared}


def rare_disease_heatmap(result: Payload, *, top: int = 30) -> dict[str, Any]:
    """Role: Provide the rare disease heatmap operation for the transforms workflow."""
    if not 1 <= top <= 30:
        raise ValueError("top must be between 1 and 30")
    records = [
        record
        for record in _records(result, "mutation_plot_data")
        if record.get("in_rare_disease")
    ]

    def disease_counts(record: Payload) -> dict[str, int]:
        """Role: Provide the disease counts operation for the transforms workflow."""
        pairs = _pairs(record.get("disease_variant_counts"))
        if pairs:
            return {
                str(disease): int(count)
                for disease, count in pairs
                if str(disease) and int(count) > 0
            }
        return {disease: 1 for disease in _strings(record.get("diseases"))}

    def rank(record: Payload) -> tuple[int, int, int, str]:
        """Role: Provide the rank operation for the transforms workflow."""
        counts = disease_counts(record)
        return (
            -sum(counts.values()),
            -len(counts),
            -(_integer(record.get("cancer_patient_count")) or 0),
            str(record.get("protein_change", "")),
        )

    shared = sorted((record for record in records if record.get("in_oncology")), key=rank)
    rare_only = sorted((record for record in records if not record.get("in_oncology")), key=rank)
    selected = (shared + rare_only)[:top]
    counts_by_mutation = {
        str(record.get("protein_change", "")): disease_counts(record) for record in selected
    }
    diseases = sorted(
        {disease for counts in counts_by_mutation.values() for disease in counts},
        key=lambda disease: (
            -sum(counts.get(disease, 0) for counts in counts_by_mutation.values()),
            disease,
        ),
    )
    mutations = [str(record.get("protein_change", "")) for record in selected]
    return {
        "mutations": mutations,
        "diseases": diseases,
        "values": [
            [counts_by_mutation[mutation].get(disease, 0) for mutation in mutations]
            for disease in diseases
        ],
        "shared": [bool(record.get("in_oncology")) for record in selected],
    }


def alphamissense_rows(result: Payload) -> list[dict[str, Any]]:
    """Role: Provide the alphamissense rows operation for the transforms workflow."""
    shared = {
        str(record.get("protein_change"))
        for record in _records(result, "shared_mutations")
    }
    rows: list[dict[str, Any]] = []
    seen: set[tuple[Any, ...]] = set()
    for record in _records(result, "evidence_records"):
        if str(record.get("source", "")).casefold() != "alphamissense":
            continue
        protein_change = str(record.get("protein_change") or "")
        score = _float(record.get("prediction_score"))
        key = (protein_change, score, record.get("prediction_class"))
        if key in seen:
            continue
        seen.add(key)
        rows.append(
            {
                "protein_change": protein_change,
                "position": _protein_position(protein_change),
                "score": score,
                "prediction": record.get("prediction_class") or "Unclassified",
                "shared_with_oncology": protein_change in shared,
                "source_record_id": record.get("source_record_id"),
                "source_version": record.get("source_version"),
            }
        )
    return sorted(
        rows,
        key=lambda row: (
            row["position"] is None,
            row["position"] or 0,
            str(row["protein_change"]),
        ),
    )


def therapy_rows(result: Payload | None) -> list[dict[str, Any]]:
    """Role: Provide the therapy rows operation for the transforms workflow."""
    if result is None:
        return []
    layer_keys = (
        ("Variant-specific evidence", "variant_therapy_evidence"),
        ("Regulatory therapy", "regulatory_therapy"),
        ("Trial and target discovery", "trial_and_target_discovery"),
    )
    rows: list[dict[str, Any]] = []
    for label, key in layer_keys:
        for record in _records(result, key):
            rows.append(
                {
                    "layer": label,
                    "therapy": record.get("therapy_name"),
                    "status": record.get("therapy_status"),
                    "protein_change": record.get("protein_change"),
                    "molecular_profile": record.get("molecular_profile"),
                    "diseases": "; ".join(_strings(record.get("disease_names"))),
                    "source": record.get("source"),
                    "evidence_level": record.get("evidence_level"),
                    "effect": record.get("effect"),
                    "trial_id": record.get("trial_id"),
                    "trial_phase": record.get("trial_phase"),
                    "jurisdiction": record.get("jurisdiction"),
                    "source_url": record.get("source_url"),
                }
            )
    return rows


def source_status_rows(*results: Payload | None) -> list[dict[str, Any]]:
    """Role: Provide the source status rows operation for the transforms workflow."""
    rows: list[dict[str, Any]] = []
    for result in results:
        if result is None:
            continue
        for record in _records(result, "source_statuses"):
            rows.append(
                {
                    "source": record.get("source"),
                    "status": record.get("status"),
                    "layer": record.get("evidence_layer"),
                    "records": _integer(record.get("record_count")) or 0,
                    "version": record.get("source_version"),
                    "message": record.get("message"),
                    "url": record.get("url"),
                }
            )
    return rows


_STRING_CHANNELS = (
    "experimental_score",
    "database_score",
    "coexpression_score",
    "textmining_score",
    "neighborhood_score",
    "fusion_score",
    "phylogenetic_score",
)


def string_partner_rows(result: Payload | None) -> list[dict[str, Any]]:
    """Return direct partners of the mapped query, retaining the strongest edge."""

    if result is None:
        return []
    query_id = str(result.get("mapped_string_id") or "")
    strongest: dict[str, dict[str, Any]] = {}
    for edge in _records(result, "edges"):
        node_a = str(edge.get("string_id_a") or "")
        node_b = str(edge.get("string_id_b") or "")
        if node_a == query_id:
            partner_id = node_b
            partner = str(edge.get("preferred_name_b") or partner_id)
        elif node_b == query_id:
            partner_id = node_a
            partner = str(edge.get("preferred_name_a") or partner_id)
        else:
            continue
        score = _float(edge.get("score")) or 0.0
        row: dict[str, Any] = {
            "partner": partner,
            "partner_string_id": partner_id,
            "combined_score": score,
        }
        for channel in _STRING_CHANNELS:
            row[channel] = _float(edge.get(channel)) or 0.0
        previous = strongest.get(partner_id)
        if previous is None or score > float(previous["combined_score"]):
            strongest[partner_id] = row
    return sorted(
        strongest.values(),
        key=lambda row: (-float(row["combined_score"]), str(row["partner"])),
    )


def string_graph_data(
    result: Payload | None,
    *,
    max_secondary_edges: int | None = None,
) -> dict[str, Any]:
    """Build a bounded graph while preserving every direct query-partner edge."""

    if result is None:
        return {"nodes": [], "edges": [], "query_id": "", "total_edges": 0}
    query_id = str(result.get("mapped_string_id") or "")
    nodes = [dict(node) for node in _records(result, "nodes")]
    edges = [dict(edge) for edge in _records(result, "edges")]
    direct = [
        edge
        for edge in edges
        if query_id
        and query_id in {str(edge.get("string_id_a")), str(edge.get("string_id_b"))}
    ]
    secondary = [edge for edge in edges if edge not in direct]
    secondary.sort(
        key=lambda edge: -(_float(edge.get("score")) or 0.0)
    )
    if max_secondary_edges is None:
        max_secondary_edges = max(8, max(0, len(nodes) - 1) // 2)
    if max_secondary_edges < 0:
        raise ValueError("max_secondary_edges must be non-negative")
    displayed = direct + secondary[:max_secondary_edges]
    degrees: dict[str, int] = defaultdict(int)
    for edge in edges:
        node_a = str(edge.get("string_id_a") or "")
        node_b = str(edge.get("string_id_b") or "")
        if node_a and node_b:
            degrees[node_a] += 1
            degrees[node_b] += 1
    for node in nodes:
        node_id = str(node.get("string_id") or "")
        node["degree"] = degrees.get(node_id, 0)
        node["is_query"] = bool(node.get("is_query")) or node_id == query_id
    return {
        "nodes": nodes,
        "edges": displayed,
        "query_id": query_id,
        "total_edges": len(edges),
        "displayed_edges": len(displayed),
    }


def structure_variant_rows(result: Payload) -> list[dict[str, Any]]:
    """List unique simple missense variants eligible for sequence modelling."""

    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for record in _records(result, "mutation_plot_data"):
        change = str(record.get("protein_change") or "").removeprefix("p.")
        match = _SIMPLE_MISSENSE.fullmatch(change)
        if change in seen or match is None:
            continue
        seen.add(change)
        rows.append(
            {
                "protein_change": change,
                "position": int(match.group(1)),
                "category": mutation_category(record),
                "cancer_patient_count": _integer(
                    record.get("cancer_patient_count")
                )
                or 0,
                "diseases": "; ".join(_strings(record.get("diseases"))),
            }
        )
    priority = {"Shared": 0, "Rare disease only": 1, "Oncology only": 2}
    return sorted(
        rows,
        key=lambda row: (
            priority[str(row["category"])],
            -int(row["cancer_patient_count"]),
            int(row["position"]),
        ),
    )


def _records(payload: Payload, key: str) -> list[Payload]:
    """Role: Compute records for the transforms workflow."""
    value = payload.get(key, [])
    if not isinstance(value, list):
        return []
    return [record for record in value if isinstance(record, dict)]


def _strings(value: Any) -> list[str]:
    """Role: Compute strings for the transforms workflow."""
    if not isinstance(value, (list, tuple)):
        return []
    return [str(item) for item in value if item is not None and str(item)]


def _pairs(value: Any) -> list[tuple[Any, Any]]:
    """Role: Compute pairs for the transforms workflow."""
    if not isinstance(value, (list, tuple)):
        return []
    pairs: list[tuple[Any, Any]] = []
    for item in value:
        if isinstance(item, (list, tuple)) and len(item) == 2:
            pairs.append((item[0], item[1]))
    return pairs


def _integer(value: Any) -> int | None:
    """Role: Compute integer for the transforms workflow."""
    if value is None or isinstance(value, bool):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _float(value: Any) -> float | None:
    """Role: Compute float for the transforms workflow."""
    if value is None or isinstance(value, bool):
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _protein_position(change: str) -> int | None:
    """Role: Compute protein position for the transforms workflow."""
    match = _POSITION.search(change)
    return int(match.group(1)) if match else None
