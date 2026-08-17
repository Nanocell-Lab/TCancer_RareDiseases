"""Command-line entry point."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from .analysis import DEFAULT_ALLOWED_CLASSIFICATIONS, find_shared_mutations
from .clients import NcbiClinVarClient
from .export import write_csv, write_json
from .plots import PlotDependencyError, generate_gene_plots


def _positive_integer(value: str) -> int:
    """Role: Parse and enforce a positive numeric command-line value."""
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("debe ser un entero positivo")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    """Role: Build the command-line parser and define supported options."""
    parser = argparse.ArgumentParser(
        description=(
            "Find ClinVar variants with ORPHAcodes and detect candidate "
            "protein matches in GDC/TCGA somatic mutations."
        )
    )
    parser.add_argument("gene", help="HGNC symbol, for example EGFR")
    parser.add_argument(
        "--project",
        action="append",
        default=[],
        help="Proyecto GDC concreto; se puede repetir (ej. TCGA-LUAD)",
    )
    parser.add_argument(
        "--program",
        default="TCGA",
        help="GDC program (default: TCGA; use ALL to disable filtering)",
    )
    parser.add_argument("--min-patients", type=int, default=1)
    parser.add_argument(
        "--exclude-condition",
        action="append",
        default=[],
        help="Nombre ClinVar a excluir; se puede repetir",
    )
    parser.add_argument(
        "--classification",
        action="append",
        default=[],
        help=(
            "Clasificacion germinal ClinVar permitida; se puede repetir. "
            "Use ALL para incluir VUS, benignas y conflictivas."
        ),
    )
    parser.add_argument(
        "--max-clinvar-records",
        type=int,
        help="Limite para smoke tests; omitir en el analisis completo",
    )
    parser.add_argument("--email", default=os.getenv("NCBI_EMAIL"))
    parser.add_argument("--api-key", default=os.getenv("NCBI_API_KEY"))
    parser.add_argument("--output-dir", default="results")
    parser.add_argument(
        "--plots",
        action="store_true",
        help="Generate lollipop, Venn, recurrence, heatmap, and dashboard plots",
    )
    parser.add_argument(
        "--plot-format",
        action="append",
        choices=("png", "svg"),
        default=[],
        help="Plot format; may be repeated (default: png and svg)",
    )
    parser.add_argument(
        "--max-plot-mutations",
        type=_positive_integer,
        default=40,
        help="Maximum mutations in lollipop and heatmap plots",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Role: Run the module entry point and coordinate its user-facing workflow."""
    args = build_parser().parse_args(argv)
    clinvar = NcbiClinVarClient(email=args.email, api_key=args.api_key)
    allowed_classifications = (
        None
        if any(value.upper() == "ALL" for value in args.classification)
        else tuple(args.classification) or DEFAULT_ALLOWED_CLASSIFICATIONS
    )
    result = find_shared_mutations(
        args.gene,
        min_cancer_patients=args.min_patients,
        program=None if args.program.upper() == "ALL" else args.program,
        projects=args.project,
        max_clinvar_records=args.max_clinvar_records,
        exclude_conditions=args.exclude_condition,
        allowed_germline_classifications=allowed_classifications,
        clinvar_client=clinvar,
    )

    output_dir = Path(args.output_dir)
    stem = result.gene.lower()
    json_path = write_json(result, output_dir / f"{stem}_report.json")
    csv_path = write_csv(result, output_dir / f"{stem}_shared.csv")
    print(
        f"{result.gene}: {result.rare_variant_count} rare variants, "
        f"{result.cancer_observation_count} oncology observations, "
        f"{result.shared_mutation_count} coincidencias candidatas."
    )
    print(f"JSON: {json_path}")
    print(f"CSV:  {csv_path}")
    if args.plots:
        try:
            plot_paths = generate_gene_plots(
                result,
                output_dir / "plots",
                formats=tuple(args.plot_format) or ("png", "svg"),
                max_lollipop_mutations=args.max_plot_mutations,
                max_heatmap_mutations=args.max_plot_mutations,
            )
        except PlotDependencyError as exc:
            raise SystemExit(str(exc)) from exc
        for name, paths in plot_paths.items():
            print(f"{name.capitalize()}: {', '.join(str(path) for path in paths)}")
    return 0
