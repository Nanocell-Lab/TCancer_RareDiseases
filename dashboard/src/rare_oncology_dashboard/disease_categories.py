"""Inherited-condition categories used by current and legacy snapshots."""

from __future__ import annotations

import re
import unicodedata
from collections.abc import Iterable, Mapping
from typing import Any


HEREDITARY_CANCER = "Hereditary cancer"
OTHER_INHERITED_DISEASE = "Other inherited diseases"
INHERITED_DISEASE_CATEGORIES = (HEREDITARY_CANCER, OTHER_INHERITED_DISEASE)

_INHERITED_CONTEXT = re.compile(
    r"\b(?:hereditar\w*|inherited|familial|germline|predispos\w*|susceptib\w*)\b"
)
_CANCER_CONTEXT = re.compile(
    r"\b(?:cancer|carcinoma|melanoma|neoplas\w*|tumou?r|sarcoma|leuk\w*|lymphoma|glioma)\b"
)
_NAMED_CANCER_SYNDROMES = (
    "li fraumeni syndrome",
    "lynch syndrome",
    "von hippel lindau syndrome",
    "cowden syndrome",
    "peutz jeghers syndrome",
    "familial adenomatous polyposis",
    "multiple endocrine neoplasia",
)


def inherited_disease_category(condition: str) -> str:
    """Classify only explicitly supported hereditary-cancer condition names."""

    normalized = _normalized_text(condition)
    if any(name in normalized for name in _NAMED_CANCER_SYNDROMES):
        return HEREDITARY_CANCER
    if _INHERITED_CONTEXT.search(normalized) and _CANCER_CONTEXT.search(
        normalized
    ):
        return HEREDITARY_CANCER
    return OTHER_INHERITED_DISEASE


def partition_inherited_diseases(
    conditions: Iterable[str],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Return unique hereditary-cancer and other inherited conditions."""

    hereditary_cancer: set[str] = set()
    other_inherited: set[str] = set()
    for condition in conditions:
        value = str(condition).strip()
        if not value:
            continue
        target = (
            hereditary_cancer
            if inherited_disease_category(value) == HEREDITARY_CANCER
            else other_inherited
        )
        target.add(value)
    return tuple(sorted(hereditary_cancer)), tuple(sorted(other_inherited))


def snapshot_disease_partition(
    record: Mapping[str, Any],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Read explicit 2.1 categories, falling back for older stored snapshots."""

    hereditary_cancer = set(_strings(record.get("hereditary_cancer_diseases")))
    other_inherited = set(_strings(record.get("other_inherited_diseases")))
    unassigned = set(_strings(record.get("diseases"))) - (
        hereditary_cancer | other_inherited
    )
    fallback_cancer, fallback_other = partition_inherited_diseases(unassigned)
    hereditary_cancer.update(fallback_cancer)
    other_inherited.update(fallback_other)
    return tuple(sorted(hereditary_cancer)), tuple(sorted(other_inherited))


def _normalized_text(value: str) -> str:
    decomposed = unicodedata.normalize("NFKD", str(value).casefold())
    ascii_text = "".join(
        character
        for character in decomposed
        if not unicodedata.combining(character)
    )
    return re.sub(r"[^a-z0-9]+", " ", ascii_text).strip()


def _strings(value: object) -> tuple[str, ...]:
    if not isinstance(value, (list, tuple)):
        return ()
    return tuple(sorted({str(item).strip() for item in value if str(item).strip()}))
