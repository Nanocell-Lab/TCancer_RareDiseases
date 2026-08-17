"""Focused CLI for generating structure-prediction inputs from a stored report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

from .storage import analysis_result_from_dict
from .structure_inputs import (
    StructureInputError,
    generate_variant_structure_inputs,
    read_protein_sequence,
)


def _positive_integer(value: str) -> int:
    """Role: Parse and enforce a positive numeric command-line value."""
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return parsed


def _nonnegative_integer(value: str) -> int:
    """Role: Parse and enforce a non-negative numeric command-line value."""
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be zero or positive")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    """Role: Build the command-line parser and define supported options."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate AlphaFold 3 JSON and Boltz-2 YAML inputs from an existing "
            "stored 2.1 molecular report."
        )
    )
    parser.add_argument("--report-json", type=Path, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--variant", action="append", default=[])
    parser.add_argument("--model-seed", action="append", type=_nonnegative_integer)
    parser.add_argument("--max-variants", type=_positive_integer, default=20)
    parser.add_argument("--protein-sequence-file", type=Path)
    parser.add_argument("--offline", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Role: Run the module entry point and coordinate its user-facing workflow."""
    args = build_parser().parse_args(argv)
    try:
        payload = _read_report(args.report_json)
        result = analysis_result_from_dict(payload)
        sequence = (
            read_protein_sequence(args.protein_sequence_file)
            if args.protein_sequence_file is not None
            else None
        )
        generated = generate_variant_structure_inputs(
            result,
            args.output_dir,
            variants=tuple(args.variant) or None,
            protein_sequence=sequence,
            sequence_database=args.database,
            offline=args.offline,
            model_seeds=tuple(args.model_seed or (1,)),
            max_variants=args.max_variants,
        )
    except (OSError, json.JSONDecodeError, TypeError, StructureInputError) as exc:
        raise SystemExit(str(exc)) from exc
    print(f"Manifest: {generated.manifest_path}")
    print(
        f"{generated.gene}: {len(generated.artifacts)} files; "
        f"{len(generated.skipped)} variants skipped."
    )
    return 0


def _read_report(path: Path) -> Mapping[str, Any]:
    """Role: Compute read report for the structure_cli workflow."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError("The molecular report must be a JSON object")
    return payload


if __name__ == "__main__":
    raise SystemExit(main())
