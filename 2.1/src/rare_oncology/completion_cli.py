"""CLI for completing missing source layers in existing 2.1 SQLite files."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from .completion import (
    DEFAULT_COMPLETION_SOURCES,
    complete_missing_queries,
    discover_completion_databases,
    write_completion_report,
)


def _positive(value: str) -> int:
    """Role: Parse and enforce a positive numeric command-line value."""
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def _bounded(minimum: int, maximum: int):
    """Role: Parse and enforce a bounded numeric command-line value."""
    def parse(value: str) -> int:
        """Role: Parse and validate a typed command-line value."""
        parsed = int(value)
        if not minimum <= parsed <= maximum:
            raise argparse.ArgumentTypeError(
                f"must be between {minimum} y {maximum}"
            )
        return parsed

    return parse


def _source_file(value: str) -> tuple[str, Path]:
    """Role: Compute source file for the completion_cli workflow."""
    source, separator, path = value.partition("=")
    if not separator or not source.strip() or not path.strip():
        raise argparse.ArgumentTypeError("use SOURCE=PATH")
    return source.strip().lower(), Path(path).expanduser()


def build_parser() -> argparse.ArgumentParser:
    """Role: Build the command-line parser and define supported options."""
    parser = argparse.ArgumentParser(
        description=(
            "Complete only missing queries for genes already stored in SQLite "
            "2.1; do not repeat ClinVar/GDC."
        )
    )
    parser.add_argument("--results-root", type=Path, default=Path("results"))
    parser.add_argument("--database", type=Path, action="append", default=[])
    parser.add_argument(
        "--source",
        action="append",
        choices=(*DEFAULT_COMPLETION_SOURCES, "oncokb"),
        default=[],
        help="Source to complete; repeat as needed. Defaults to all open sources.",
    )
    parser.add_argument(
        "--regulatory-file", action="append", type=_source_file, default=[],
        metavar="SOURCE=PATH",
    )
    parser.add_argument("--oncokb-token", default=os.getenv("ONCOKB_TOKEN"))
    parser.add_argument("--max-therapy-records", type=_positive, default=100)
    parser.add_argument("--gnomad-max-variants", type=_bounded(1, 25), default=25)
    parser.add_argument("--string-score", type=_bounded(0, 1000), default=700)
    parser.add_argument("--string-partners", type=_bounded(0, 100), default=20)
    parser.add_argument(
        "--string-network-type", choices=("functional", "physical"),
        default="functional",
    )
    parser.add_argument("--max-genes", type=_positive)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--retry-errors", action="store_true")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/completion")
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Role: Run the module entry point and coordinate its user-facing workflow."""
    args = build_parser().parse_args(argv)
    databases = tuple(args.database) or discover_completion_databases(
        args.results_root
    )
    if not databases:
        raise SystemExit("No SQLite 2.1 databases with execution_cache were found")
    selected = tuple(args.source) or DEFAULT_COMPLETION_SOURCES
    try:
        report = complete_missing_queries(
            databases,
            sources=selected,
            oncokb_token=args.oncokb_token,
            regulatory_files=args.regulatory_file,
            max_therapy_records=args.max_therapy_records,
            gnomad_max_variants=args.gnomad_max_variants,
            string_score=args.string_score,
            string_partners=args.string_partners,
            string_network_type=args.string_network_type,
            offline=args.offline,
            dry_run=args.dry_run,
            retry_errors=args.retry_errors,
            strict=args.strict,
            max_genes=args.max_genes,
        )
        paths = write_completion_report(report, args.output_dir)
    except (OSError, TypeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    print(
        f"Completeness: {report.gene_count} gene/SQLite entries; "
        f"{report.query_count} queries executed, "
        f"{report.cache_hit_count} reused, "
        f"{report.missing_count} pending and {report.error_count} errors."
    )
    print("Report: " + ", ".join(str(path) for path in paths))
    return 1 if report.error_count and args.strict else 0


if __name__ == "__main__":
    raise SystemExit(main())
