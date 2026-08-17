"""CLI for incremental cross-gene protein-family analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

from .global_analysis import (
    build_global_family_analysis,
    discover_global_databases,
    write_global_report,
)
from .global_plots import generate_global_family_plots
from .plots import PlotDependencyError


def build_parser() -> argparse.ArgumentParser:
    """Role: Build the command-line parser and define supported options."""
    parser = argparse.ArgumentParser(
        description=(
            "Incrementally aggregate SQLite 2.1 snapshots and compare actionable "
            "protein families in oncology without querying the web."
        )
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=Path("results"),
        help="Root used to discover SQLite databases when --database is omitted.",
    )
    parser.add_argument(
        "--database",
        action="append",
        type=Path,
        default=[],
        help="Source SQLite database; may be repeated.",
    )
    parser.add_argument(
        "--index-database",
        type=Path,
        default=Path("results/global_family_index.sqlite3"),
        help="Derived incremental index (does not replace source databases).",
    )
    parser.add_argument(
        "--gene-list-dir",
        type=Path,
        default=Path("gene_lists"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/global_family"),
    )
    parser.add_argument("--plots", action="store_true")
    parser.add_argument(
        "--plot-format",
        action="append",
        choices=("png", "svg", "pdf"),
        default=[],
    )
    parser.add_argument("--plot-dpi", type=int, default=200)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Role: Run the module entry point and coordinate its user-facing workflow."""
    args = build_parser().parse_args(argv)
    databases = tuple(args.database) or discover_global_databases(
        args.results_root, exclude=(args.index_database,)
    )
    if not databases:
        raise SystemExit("No SQLite 2.1 databases were found for aggregation")
    try:
        report = build_global_family_analysis(
            databases,
            index_database=args.index_database,
            gene_list_dir=args.gene_list_dir,
        )
        table_paths = write_global_report(report, args.output_dir)
        plot_paths: tuple[Path, ...] = ()
        if args.plots:
            plot_paths = generate_global_family_plots(
                report,
                args.output_dir / "plots",
                formats=tuple(args.plot_format) or ("png", "svg"),
                dpi=args.plot_dpi,
            )
    except (OSError, ValueError, PlotDependencyError) as exc:
        raise SystemExit(str(exc)) from exc

    print(
        f"Global analysis: {report.active_analysis_count} genes, "
        f"{len(report.family_summaries)} subfamilies, "
        f"{report.indexed_new_snapshots} new snapshots indexed."
    )
    print(
        f"Compatible therapy snapshots: {report.active_therapy_count}; "
        f"genes with an incompatible therapy snapshot: "
        f"{report.unmatched_therapy_gene_count}."
    )
    if report.index_error_count or report.source_error_count:
        print(
            f"Warning: {report.index_error_count} snapshots skipped and "
            f"{report.source_error_count} databases with read errors."
        )
    print("Tables: " + ", ".join(str(path) for path in table_paths))
    if plot_paths:
        print("Figures: " + ", ".join(str(path) for path in plot_paths))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
