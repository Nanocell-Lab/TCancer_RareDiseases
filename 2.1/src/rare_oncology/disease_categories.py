"""Deterministic condition-name categories for inherited-disease results."""

from __future__ import annotations

import re
import unicodedata
from collections.abc import Iterable


HEREDITARY_CANCER = "Hereditary cancer"
OTHER_INHERITED_DISEASE = "Other inherited disease"

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
    """Classify a ClinVar condition without conflating it with somatic oncology.

    ClinVar does not expose a stable hereditary-cancer category in the stored
    snapshots, so classification is deliberately based on explicit wording in
    the condition name. Ambiguous names remain in the general inherited-disease
    category rather than being promoted to hereditary cancer.
    """

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


def _normalized_text(value: str) -> str:
    decomposed = unicodedata.normalize("NFKD", str(value).casefold())
    ascii_text = "".join(
        character
        for character in decomposed
        if not unicodedata.combining(character)
    )
    return re.sub(r"[^a-z0-9]+", " ", ascii_text).strip()
