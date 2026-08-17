"""Focused CLI for fetching one human STRING network without rerunning analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

from .plots import PlotDependencyError, generate_string_network_plots
from .string_network import (
    StringNetworkError,
    get_string_network,
    write_string_edges_csv,
    write_string_network_json,
    write_string_nodes_csv,
)


def _bounded_integer(minimum: int, maximum: int):
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


def build_parser() -> argparse.ArgumentParser:
    """Role: Build the command-line parser and define supported options."""
    parser = argparse.ArgumentParser(
        description="Retrieve a human STRING v12.0 network with SQLite caching."
    )
    parser.add_argument("gene", help="HGNC symbol, for example EGFR")
    parser.add_argument("--database", type=Path, default=Path("results/rare_oncology.sqlite3"))
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument(
        "--string-score", type=_bounded_integer(0, 1000), default=700
    )
    parser.add_argument(
        "--string-partners", type=_bounded_integer(0, 100), default=20
    )
    parser.add_argument(
        "--string-network-type",
        choices=("functional", "physical"),
        default="functional",
    )
    parser.add_argument("--plots", action="store_true")
    parser.add_argument(
        "--plot-format",
        action="append",
        choices=("png", "svg"),
        default=[],
    )
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--offline", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Role: Run the module entry point and coordinate its user-facing workflow."""
    args = build_parser().parse_args(argv)
    if args.refresh and args.offline:
        raise SystemExit("--refresh and --offline cannot be combined")
    try:
        run = get_string_network(
            args.gene,
            database=args.database,
            required_score=args.string_score,
            network_type=args.string_network_type,
            added_nodes=args.string_partners,
            refresh=args.refresh,
            offline=args.offline,
        )
    except (StringNetworkError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    network = run.network
    stem = network.gene.lower()
    json_path = write_string_network_json(
        network, args.output_dir / f"{stem}_string_network.json"
    )
    nodes_path = write_string_nodes_csv(
        network, args.output_dir / f"{stem}_string_nodes.csv"
    )
    edges_path = write_string_edges_csv(
        network, args.output_dir / f"{stem}_string_edges.csv"
    )
    print(
        f"{network.gene}: STRING v{network.string_version}, "
        f"{len(network.nodes)} nodes, {len(network.edges)} edges "
        f"[{'HIT SQLite' if run.cache_hit else 'WEB→SQLite'}]."
    )
    print(f"JSON: {json_path}")
    print(f"Nodes: {nodes_path}")
    print(f"Edges: {edges_path}")
    if args.plots:
        try:
            paths = generate_string_network_plots(
                network,
                args.output_dir / "plots",
                formats=tuple(args.plot_format) or ("png", "svg"),
            )
        except PlotDependencyError as exc:
            raise SystemExit(str(exc)) from exc
        print("Plots: " + ", ".join(str(path) for path in paths))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
