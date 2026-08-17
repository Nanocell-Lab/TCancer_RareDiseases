"""Build normalized local FDA and EMA tables for regulatory enrichment.

The source parsers are intentionally conservative:

* FDA pharmacogenomic-table entries are labelled ``label_biomarker`` because
  presence in that table does not by itself prove biomarker-selected approval.
* EMA entries are emitted only when the therapeutic indication explicitly
  mentions a requested gene symbol or a small, controlled biomarker alias.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import shlex
import sys
import tempfile
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, Iterable, Sequence


FDA_URL = (
    "https://www.fda.gov/drugs/science-and-research-drugs/"
    "table-pharmacogenomic-biomarkers-drug-labeling/"
)
EMA_URL = (
    "https://www.ema.europa.eu/en/documents/report/"
    "medicines-output-medicines_json-report_en.json"
)

FIELDNAMES = (
    "gene",
    "source_record_id",
    "therapy_name",
    "therapy_id",
    "protein_change",
    "biomarker_text",
    "disease_id",
    "disease_name",
    "therapy_status",
    "effect",
    "evidence_type",
    "evidence_level",
    "jurisdiction",
    "citation_id",
    "source_url",
    "source_version",
    "license",
    "therapeutic_area",
    "active_substance",
    "regulatory_source_status",
    "description",
    "retrieved_at",
)

# These aliases occur in regulatory indications and are sufficiently specific
# for candidate biomarker matching. They do not change the source's assertion.
GENE_ALIASES: dict[str, tuple[str, ...]] = {
    "ERBB2": (r"HER[ -]?2",),
    "MET": (
        r"METex14",
        r"mesenchymal[ -]epithelial transition factor gene exon 14",
    ),
    "NTRK1": (r"NTRK gene fusion",),
    "NTRK2": (r"NTRK gene fusion",),
    "NTRK3": (r"NTRK gene fusion",),
}


@dataclass(frozen=True)
class FdaRow:
    drug: str
    therapeutic_area: str
    biomarker: str
    labeling_sections: str
    source_url: str


class _FdaTableParser(HTMLParser):
    def __init__(self, base_url: str) -> None:
        """Role: Compute init for the regulatory_tables workflow."""
        super().__init__(convert_charrefs=True)
        self.base_url = base_url
        self.rows: list[FdaRow] = []
        self._in_table = False
        self._in_row = False
        self._in_cell = False
        self._cells: list[list[str]] = []
        self._hrefs: list[str | None] = []

    def handle_starttag(
        self, tag: str, attrs: list[tuple[str, str | None]]
    ) -> None:
        """Role: Provide the handle starttag operation for the regulatory_tables workflow."""
        attributes = dict(attrs)
        if tag == "table" and "pharmacogenomic" in (
            attributes.get("summary") or ""
        ).casefold():
            self._in_table = True
            return
        if not self._in_table:
            return
        if tag == "tr":
            self._in_row = True
            self._cells = []
            self._hrefs = []
        elif tag in {"td", "th"} and self._in_row:
            self._in_cell = True
            self._cells.append([])
            self._hrefs.append(None)
        elif tag == "br" and self._in_cell and self._cells:
            self._cells[-1].append(" ")
        elif tag == "a" and self._in_cell and self._hrefs:
            href = attributes.get("href")
            if href and self._hrefs[-1] is None:
                self._hrefs[-1] = urllib.parse.urljoin(self.base_url, href)

    def handle_endtag(self, tag: str) -> None:
        """Role: Provide the handle endtag operation for the regulatory_tables workflow."""
        if not self._in_table:
            return
        if tag in {"td", "th"}:
            self._in_cell = False
        elif tag == "tr" and self._in_row:
            self._finish_row()
            self._in_row = False
        elif tag == "table":
            self._in_table = False

    def handle_data(self, data: str) -> None:
        """Role: Provide the handle data operation for the regulatory_tables workflow."""
        if self._in_cell and self._cells:
            self._cells[-1].append(data)

    def _finish_row(self) -> None:
        """Role: Compute finish row for the regulatory_tables workflow."""
        cells = [_clean_text("".join(parts)) for parts in self._cells]
        if len(cells) != 4 or cells[0].casefold() == "drug":
            return
        source_url = self._hrefs[0] or self.base_url
        if source_url.startswith("http://www.accessdata.fda.gov/"):
            source_url = "https://" + source_url.removeprefix("http://")
        self.rows.append(
            FdaRow(
                drug=cells[0],
                therapeutic_area=cells[1],
                biomarker=cells[2],
                labeling_sections=cells[3],
                source_url=source_url,
            )
        )


def parse_fda_html(
    content: str,
    genes: Iterable[str],
    *,
    source_url: str = FDA_URL,
    retrieved_at: str | None = None,
) -> list[dict[str, str]]:
    """Parse the official FDA pharmacogenomic-label table."""

    gene_symbols = _normalize_genes(genes)
    parser = _FdaTableParser(source_url)
    parser.feed(content)
    if not parser.rows:
        raise ValueError("The FDA pharmacogenomic table was not found")
    version = _fda_source_version(content)
    retrieval = retrieved_at or _utc_now()
    output: list[dict[str, str]] = []
    for row_number, row in enumerate(parser.rows, start=1):
        matches = _matched_genes(row.biomarker, gene_symbols, aliases=False)
        for gene in matches:
            application = _fda_application_number(row.source_url)
            record_suffix = application or str(row_number)
            therapy_name = re.sub(r"\s+\(\d+\)$", "", row.drug).strip()
            output.append(
                _row(
                    gene=gene,
                    source_record_id=f"FDA-PGX:{record_suffix}:{gene}:{row_number}",
                    therapy_name=therapy_name,
                    therapy_id=application,
                    biomarker_text=row.biomarker,
                    therapy_status="label_biomarker",
                    evidence_type=row.labeling_sections,
                    jurisdiction="US",
                    source_url=row.source_url,
                    source_version=version,
                    therapeutic_area=row.therapeutic_area,
                    regulatory_source_status="pharmacogenomic_label_table",
                    description=(
                        f"FDA pharmacogenomic labeling entry; therapeutic area: "
                        f"{row.therapeutic_area}"
                    ),
                    retrieved_at=retrieval,
                )
            )
    return _sorted_unique(output)


def parse_ema_json(
    payload: dict[str, Any],
    genes: Iterable[str],
    *,
    source_url: str = EMA_URL,
    retrieved_at: str | None = None,
) -> list[dict[str, str]]:
    """Parse EMA's official medicine-pages JSON export."""

    raw_records = payload.get("data")
    if not isinstance(raw_records, list):
        raise ValueError("The EMA JSON does not contain a data list")
    meta = payload.get("meta") or {}
    version = _text(meta.get("timestamp"))
    retrieval = retrieved_at or _utc_now()
    gene_symbols = _normalize_genes(genes)
    output: list[dict[str, str]] = []
    for index, item in enumerate(raw_records, start=1):
        if not isinstance(item, dict) or _text(item.get("category")).casefold() != "human":
            continue
        indication = _clean_text(_text(item.get("therapeutic_indication")))
        matched = _matched_genes(indication, gene_symbols, aliases=True)
        matches = tuple(
            (gene, effect)
            for gene in matched
            if (effect := _ema_biomarker_effect(indication, gene)) is not None
        )
        if not matches:
            continue
        medicine = _clean_text(_text(item.get("name_of_medicine")))
        product_number = _clean_text(_text(item.get("ema_product_number")))
        medicine_url = _clean_text(_text(item.get("medicine_url"))) or source_url
        if not medicine or not product_number:
            continue
        source_status = _clean_text(_text(item.get("medicine_status")))
        status = _ema_status(source_status, item)
        source_version = version or _clean_text(_text(item.get("last_updated_date")))
        for gene, effect in matches:
            output.append(
                _row(
                    gene=gene,
                    source_record_id=f"{product_number}:{gene}",
                    therapy_name=medicine,
                    therapy_id=product_number,
                    biomarker_text=_biomarker_excerpt(indication, gene),
                    disease_name=_clean_text(
                        _text(item.get("therapeutic_area_mesh"))
                    ),
                    therapy_status=status,
                    effect=effect,
                    evidence_type="EMA therapeutic indication",
                    evidence_level=_ema_evidence_level(item),
                    jurisdiction="EU",
                    source_url=medicine_url,
                    source_version=source_version,
                    therapeutic_area=_clean_text(
                        _text(item.get("therapeutic_area_mesh"))
                    ),
                    active_substance=_clean_text(_text(item.get("active_substance"))),
                    regulatory_source_status=source_status,
                    description=indication,
                    retrieved_at=retrieval,
                )
            )
    return _sorted_unique(output)


def write_table(path: str | Path, rows: Iterable[dict[str, str]]) -> Path:
    """Atomically write a normalized, spreadsheet-safe UTF-8 CSV."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=output.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _safe_csv_value(row.get(key, "")) for key in FIELDNAMES})
    temporary.replace(output)
    return output


def load_genes(paths: Iterable[str | Path], inline: str | None = None) -> tuple[str, ...]:
    """Role: Load and validate genes data."""
    genes: set[str] = set()
    for path in paths:
        source = Path(path)
        if not source.is_file():
            raise FileNotFoundError(f"Gene list does not exist: {source}")
        for line in source.read_text(encoding="utf-8").splitlines():
            value = line.split("#", 1)[0].strip()
            if value:
                genes.add(_validate_gene(value))
    if inline:
        for value in inline.split(","):
            if value.strip():
                genes.add(_validate_gene(value.strip()))
    if not genes:
        raise ValueError("At least one gene must be provided")
    return tuple(sorted(genes))


def build_parser() -> argparse.ArgumentParser:
    """Role: Build the command-line parser and define supported options."""
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description=(
            "Download FDA/EMA data and generate normalized regulatory tables for "
            "rare-oncology 2.1."
        )
    )
    parser.add_argument(
        "--source", choices=("all", "fda", "ema"), default="all"
    )
    parser.add_argument(
        "--output-dir", default=str(root / "downloads" / "regulatory")
    )
    parser.add_argument(
        "--gene-list",
        action="append",
        default=None,
        help="HGNC symbol list; may be repeated.",
    )
    parser.add_argument("--genes", help="Additional symbols separated by commas.")
    parser.add_argument("--fda-input", help="HTML FDA local; evita descargar FDA.")
    parser.add_argument("--ema-input", help="JSON EMA local; evita descargar EMA.")
    parser.add_argument("--fda-url", default=FDA_URL)
    parser.add_argument("--ema-url", default=EMA_URL)
    parser.add_argument("--timeout", type=_positive_float, default=90.0)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Role: Run the module entry point and coordinate its user-facing workflow."""
    args = build_parser().parse_args(argv)
    root = Path(__file__).resolve().parents[2]
    default_lists = (
        root / "gene_lists" / "receptor_tyrosine_kinases.txt",
        root / "gene_lists" / "tgfb_receptors.txt",
    )
    try:
        genes = load_genes(args.gene_list or default_lists, args.genes)
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        environment: dict[str, Path] = {}
        if args.source in {"all", "fda"}:
            fda_bytes = _read_or_download(
                args.fda_input, args.fda_url, args.timeout
            )
            raw_fda = output_dir / "fda_pharmacogenomic_biomarkers.html"
            raw_fda.write_bytes(fda_bytes)
            fda_rows = parse_fda_html(
                fda_bytes.decode("utf-8", "replace"), genes, source_url=args.fda_url
            )
            fda_output = write_table(output_dir / "fda_biomarkers.csv", fda_rows)
            environment["FDA_FILE"] = fda_output.resolve()
            print(f"FDA: {len(fda_rows)} rows -> {fda_output}")
        if args.source in {"all", "ema"}:
            ema_bytes = _read_or_download(
                args.ema_input, args.ema_url, args.timeout
            )
            raw_ema = output_dir / "ema_medicines.json"
            raw_ema.write_bytes(ema_bytes)
            ema_payload = json.loads(ema_bytes)
            if not isinstance(ema_payload, dict):
                raise ValueError("The EMA JSON root is not an object")
            ema_rows = parse_ema_json(ema_payload, genes, source_url=args.ema_url)
            ema_output = write_table(output_dir / "ema_labels.csv", ema_rows)
            environment["EMA_FILE"] = ema_output.resolve()
            print(f"EMA: {len(ema_rows)} rows -> {ema_output}")
        env_path = output_dir / "regulatory_tables.env"
        env_path.write_text(
            "".join(
                f"export {name}={shlex.quote(str(path))}\n"
                for name, path in sorted(environment.items())
            ),
            encoding="utf-8",
        )
        print(f"Variables: source {env_path}")
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    return 0


def _read_or_download(local_path: str | None, url: str, timeout: float) -> bytes:
    """Role: Compute read or download for the regulatory_tables workflow."""
    if local_path:
        return Path(local_path).read_bytes()
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json,text/html;q=0.9,*/*;q=0.1",
            "User-Agent": "rare-oncology-regulatory-table-builder/2.1",
        },
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def _matched_genes(
    text: str, genes: tuple[str, ...], *, aliases: bool
) -> tuple[str, ...]:
    """Role: Compute matched genes for the regulatory_tables workflow."""
    matches: list[str] = []
    for gene in genes:
        if re.search(_gene_search_pattern(gene, aliases=aliases), text, re.IGNORECASE):
            matches.append(gene)
    return tuple(matches)


def _gene_search_pattern(gene: str, *, aliases: bool) -> str:
    """Role: Compute gene search pattern for the regulatory_tables workflow."""
    patterns = [rf"(?<![A-Za-z0-9]){re.escape(gene)}(?![A-Za-z0-9])"]
    if aliases:
        patterns.extend(GENE_ALIASES.get(gene, ()))
    return "(?:" + "|".join(patterns) + ")"


def _biomarker_excerpt(indication: str, gene: str) -> str:
    """Role: Compute biomarker excerpt for the regulatory_tables workflow."""
    sentences = re.split(r"(?<=[.!?])\s+", indication)
    for sentence in sentences:
        if _matched_genes(sentence, (gene,), aliases=True):
            return sentence[:2000]
    return indication[:2000]


def _ema_biomarker_effect(indication: str, gene: str) -> str | None:
    """Return a conservative biomarker role, excluding prior-therapy mentions."""

    for sentence in re.split(r"(?<=[.!?])\s+", indication):
        if not _matched_genes(sentence, (gene,), aliases=True):
            continue
        without_prior_therapy = re.sub(
            rf"\banti[- ]?(?:{re.escape(gene)}|HER[ -]?2|NTRK)\b",
            "",
            sentence,
            flags=re.IGNORECASE,
        )
        if not _matched_genes(without_prior_therapy, (gene,), aliases=True):
            continue
        if not re.search(
            r"\b(?:biomarker|positive|negative|mutation|mutant|alteration|"
            r"deletion|insertion|amplif\w*|fusion|rearrang\w*|aberration|"
            r"express\w*|overexpress\w*|status|wild[- ]?type|exon)\b",
            sentence,
            flags=re.IGNORECASE,
        ):
            continue
        symbol_pattern = _gene_search_pattern(gene, aliases=True)
        negative_before = re.search(
            rf"\b(?:no|without)\b[^.;:]{{0,50}}(?:{symbol_pattern})",
            sentence,
            flags=re.IGNORECASE,
        )
        negative_after = re.search(
            rf"(?:{symbol_pattern})[^.;:]{{0,35}}\bnegative\b",
            sentence,
            flags=re.IGNORECASE,
        )
        return (
            "biomarker_exclusion"
            if negative_before or negative_after
            else "biomarker_selected"
        )
    return None


def _ema_status(source_status: str, item: dict[str, Any]) -> str:
    """Role: Compute ema status for the regulatory_tables workflow."""
    status = source_status.casefold()
    if status == "authorised":
        if _text(item.get("conditional_approval")).casefold() == "yes":
            return "conditional_approval"
        return "regulatory_approved"
    mappings = {
        "withdrawn": "withdrawn",
        "application withdrawn": "application_withdrawn",
        "withdrawn from rolling review": "rolling_review_withdrawn",
        "refused": "refused",
        "opinion": "regulatory_opinion",
        "opinion under re-examination": "opinion_under_reexamination",
        "lapsed": "lapsed",
        "expired": "expired",
        "revoked": "revoked",
        "suspended": "suspended",
    }
    return mappings.get(status, "regulatory_status_not_provided")


def _ema_evidence_level(item: dict[str, Any]) -> str:
    """Role: Compute ema evidence level for the regulatory_tables workflow."""
    flags = []
    for key, label in (
        ("conditional_approval", "conditional approval"),
        ("exceptional_circumstances", "exceptional circumstances"),
        ("orphan_medicine", "orphan medicine"),
        ("prime_priority_medicine", "PRIME"),
        ("accelerated_assessment", "accelerated assessment"),
    ):
        if _text(item.get(key)).casefold() == "yes":
            flags.append(label)
    return " | ".join(flags)


def _fda_source_version(content: str) -> str:
    """Role: Compute fda source version for the regulatory_tables workflow."""
    match = re.search(
        r"Content current as of:.*?<time[^>]*datetime=[\"']([^\"']+)",
        content,
        flags=re.IGNORECASE | re.DOTALL,
    )
    return match.group(1) if match else ""


def _fda_application_number(url: str) -> str:
    """Role: Compute fda application number for the regulatory_tables workflow."""
    query = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
    values = query.get("varApplNo") or query.get("varapplno") or []
    return values[0] if values else ""


def _row(**values: str) -> dict[str, str]:
    """Role: Compute row for the regulatory_tables workflow."""
    return {field: values.get(field, "") for field in FIELDNAMES}


def _sorted_unique(rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    """Role: Compute sorted unique for the regulatory_tables workflow."""
    keyed: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in rows:
        key = (row["gene"], row["source_record_id"], row["therapy_name"])
        keyed.setdefault(key, row)
    return sorted(
        keyed.values(),
        key=lambda row: (
            row["gene"],
            row["therapy_name"].casefold(),
            row["source_record_id"],
        ),
    )


def _normalize_genes(genes: Iterable[str]) -> tuple[str, ...]:
    """Role: Normalize genes into a stable internal representation."""
    return tuple(sorted({_validate_gene(value) for value in genes}))


def _validate_gene(value: str) -> str:
    """Role: Validate gene constraints before execution."""
    gene = value.strip().upper()
    if not re.fullmatch(r"[A-Z][A-Z0-9-]*", gene):
        raise ValueError(f"Invalid gene symbol: {value!r}")
    return gene


def _clean_text(value: str) -> str:
    """Role: Compute clean text for the regulatory_tables workflow."""
    return re.sub(r"\s+", " ", html.unescape(value)).strip()


def _text(value: object) -> str:
    """Role: Compute text for the regulatory_tables workflow."""
    return "" if value is None else str(value)


def _utc_now() -> str:
    """Role: Return the current UTC timestamp in canonical form."""
    return datetime.now(timezone.utc).isoformat()


def _positive_float(value: str) -> float:
    """Role: Parse and enforce a positive numeric command-line value."""
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number") from exc
    if not 0 < parsed < float("inf"):
        raise argparse.ArgumentTypeError("must be finite and greater than zero")
    return parsed


def _safe_csv_value(value: str) -> str:
    """Role: Sanitize csv value for safe serialization."""
    trimmed = value.lstrip(" \t\r\n")
    if value.startswith(("\t", "\r", "\n")) or trimmed.startswith(
        ("=", "+", "-", "@")
    ):
        return "'" + value
    return value


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
