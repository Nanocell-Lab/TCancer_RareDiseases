"""Conservative protein-change normalization for cross-source comparison."""

from __future__ import annotations

import re


_AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
    "Stop": "*",
}

_VALID_PROTEIN_CHANGE = re.compile(
    r"^[A-Z*]\d+(?:[A-Z*=]|_[A-Z*]\d+(?:[A-Z*=])?)"
    r"(?:fs\*?\d*|del(?:ins[A-Z*]+)?|dup|ins[A-Z*]+)?$"
)


def normalize_protein_change(value: str | None) -> str | None:
    """Return a comparable one-letter protein change, or ``None`` if uncertain.

    The function deliberately rejects unknown values and does not infer transcript
    equivalence. Examples: ``p.Leu858Arg`` -> ``L858R`` and
    ``p.Glu746_Ala750del`` -> ``E746_A750del``.
    """

    if not value:
        return None

    change = value.strip().strip("()")
    if change.startswith("p."):
        change = change[2:]
    change = change.replace(" ", "")
    if not change or "?" in change:
        return None

    for three, one in _AA3_TO_1.items():
        change = re.sub(three, one, change, flags=re.IGNORECASE)

    change = change.replace("X", "*")
    if not _VALID_PROTEIN_CHANGE.fullmatch(change):
        return None
    return change


def primary_protein_change(summary: dict) -> tuple[str, str] | None:
    """Extract ClinVar's title-level protein change, with a summary fallback."""

    title = str(summary.get("title", ""))
    title_match = re.search(r"\(p\.([^)]+)\)", title)
    candidates: list[str] = []
    if title_match:
        candidates.append(title_match.group(1))
    candidates.extend(str(summary.get("protein_change", "")).split(","))

    for raw in candidates:
        normalized = normalize_protein_change(raw)
        if normalized:
            return normalized, raw.strip()
    return None


def clinvar_spdi_key(value: str | None) -> tuple[str, int, str, str] | None:
    """Convert a simple GRCh38 SPDI substitution to a genomic comparison key."""

    if not value:
        return None
    match = re.fullmatch(r"NC_0*(\d+)\.\d+:(\d+):([ACGT]+):([ACGT]+)", value)
    if not match or len(match.group(3)) != 1 or len(match.group(4)) != 1:
        return None
    chromosome_number = int(match.group(1))
    chromosome = {23: "X", 24: "Y", 12920: "MT"}.get(
        chromosome_number, str(chromosome_number)
    )
    # SPDI uses a zero-based interbase position; GDC g.HGVS is one-based.
    return (f"chr{chromosome}", int(match.group(2)) + 1, match.group(3), match.group(4))


def gdc_genomic_change_key(value: str | None) -> tuple[str, int, str, str] | None:
    """Convert a simple GDC genomic substitution to a comparison key."""

    if not value:
        return None
    match = re.fullmatch(r"(chr(?:\d+|X|Y|M|MT)):g\.(\d+)([ACGT])>([ACGT])", value)
    if not match:
        return None
    chromosome = "chrMT" if match.group(1) == "chrM" else match.group(1)
    return (chromosome, int(match.group(2)), match.group(3), match.group(4))


def protein_position(value: str | None) -> int | None:
    """Extract the first amino-acid position from a normalized protein change."""

    normalized = normalize_protein_change(value)
    if normalized is None:
        return None
    match = re.match(r"^[A-Z*](\d+)", normalized)
    return int(match.group(1)) if match else None
