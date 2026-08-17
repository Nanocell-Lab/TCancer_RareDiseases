"""Retrieve curated protein domains and secondary structure from UniProt."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from .http import JsonHttpClient


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


@dataclass(frozen=True)
class ProteinFeature:
    kind: str
    start: int
    end: int
    label: str


@dataclass(frozen=True)
class ProteinAnnotation:
    gene: str
    accession: str | None
    length: int
    domains: tuple[ProteinFeature, ...]
    secondary_structure: tuple[ProteinFeature, ...]
    source_url: str | None


class UniProtProteinClient:
    """Fetch the reviewed canonical human UniProt entry for one gene."""

    def __init__(self, http: JsonHttpClient | None = None) -> None:
        """Role: Compute init for the protein_annotations workflow."""
        self.http = http or JsonHttpClient(
            user_agent="rare-oncology/2.0 protein-annotation",
            retries=3,
        )

    def fetch(self, gene: str) -> ProteinAnnotation | None:
        """Role: Retrieve and normalize the requested remote resource."""
        normalized_gene = gene.strip().upper()
        payload = self.http.get_json(
            UNIPROT_SEARCH_URL,
            {
                "query": (
                    f"(gene_exact:{normalized_gene}) AND "
                    "(organism_id:9606) AND (reviewed:true)"
                ),
                "format": "json",
                "size": 1,
            },
        )
        results = payload.get("results")
        if not isinstance(results, list) or not results:
            return None
        entry = results[0]
        if not isinstance(entry, dict):
            return None
        accession = _string_or_none(entry.get("primaryAccession"))
        sequence = entry.get("sequence")
        length = _positive_int(sequence.get("length")) if isinstance(sequence, dict) else None
        if length is None:
            return None

        domains: list[ProteinFeature] = []
        secondary: list[ProteinFeature] = []
        features = entry.get("features")
        if isinstance(features, list):
            for feature in features:
                parsed = _parse_feature(feature)
                if parsed is None:
                    continue
                if parsed.kind == "domain":
                    domains.append(parsed)
                elif parsed.kind in {"helix", "strand"}:
                    secondary.append(parsed)

        return ProteinAnnotation(
            gene=normalized_gene,
            accession=accession,
            length=length,
            domains=tuple(sorted(domains, key=_feature_sort_key)),
            secondary_structure=tuple(sorted(secondary, key=_feature_sort_key)),
            source_url=(
                f"https://rest.uniprot.org/uniprotkb/{accession}.json"
                if accession
                else None
            ),
        )


@lru_cache(maxsize=128)
def fetch_protein_annotation(gene: str) -> ProteinAnnotation | None:
    """Role: Retrieve and normalize protein annotation data."""
    return UniProtProteinClient().fetch(gene)


def inferred_annotation(gene: str, positions: list[int]) -> ProteinAnnotation:
    """Role: Provide the inferred annotation operation for the protein_annotations workflow."""
    return ProteinAnnotation(
        gene=gene,
        accession=None,
        length=max(positions, default=1),
        domains=(),
        secondary_structure=(),
        source_url=None,
    )


def _parse_feature(value: object) -> ProteinFeature | None:
    """Role: Parse and validate feature from its source representation."""
    if not isinstance(value, dict):
        return None
    feature_type = value.get("type")
    kind_by_type = {
        "Domain": "domain",
        "Helix": "helix",
        "Beta strand": "strand",
    }
    kind = kind_by_type.get(feature_type) if isinstance(feature_type, str) else None
    location = value.get("location")
    if kind is None or not isinstance(location, dict):
        return None
    start_value = location.get("start")
    end_value = location.get("end")
    start = _location_value(start_value)
    end = _location_value(end_value)
    if start is None or end is None or end < start:
        return None
    description = _string_or_none(value.get("description"))
    labels = {
        "domain": description or "Domain",
        "helix": "Helix",
        "strand": "β-strand",
    }
    return ProteinFeature(kind=kind, start=start, end=end, label=labels[kind])


def _location_value(value: object) -> int | None:
    """Role: Compute location value for the protein_annotations workflow."""
    if not isinstance(value, dict):
        return None
    return _positive_int(value.get("value"))


def _positive_int(value: object) -> int | None:
    """Role: Parse and enforce a positive numeric command-line value."""
    if isinstance(value, bool):
        return None
    try:
        parsed = int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return parsed if parsed > 0 else None


def _string_or_none(value: object) -> str | None:
    """Role: Compute string or none for the protein_annotations workflow."""
    return value.strip() if isinstance(value, str) and value.strip() else None


def _feature_sort_key(feature: ProteinFeature) -> tuple[int, int, str]:
    """Role: Compute feature sort key for the protein_annotations workflow."""
    return feature.start, feature.end, feature.label
