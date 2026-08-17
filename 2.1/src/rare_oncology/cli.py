"""Command-line entry point."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

from .alphamissense import AlphaMissenseEvidenceProvider
from .analysis import DEFAULT_ALLOWED_CLASSIFICATIONS, find_shared_mutations
from .clients import NcbiClinVarClient
from .evidence import (
    GnomadEvidenceProvider,
    LocalTableEvidenceProvider,
    OrphadataXmlEvidenceProvider,
    PanelAppEvidenceProvider,
)
from .export import (
    write_csv,
    write_evidence_csv,
    write_json,
    write_therapy_csv,
    write_therapy_json,
)
from .plots import (
    PlotDependencyError,
    generate_gene_plots,
    generate_string_network_plots,
)
from .protein_annotations import inferred_annotation
from .storage import (
    CacheMissError,
    find_shared_mutations_cached,
    find_therapy_associations_cached,
)
from .structure_inputs import (
    StructureInputError,
    generate_variant_structure_inputs,
    read_protein_sequence,
)
from .string_network import (
    StringNetworkError,
    get_string_network,
    write_string_edges_csv,
    write_string_network_json,
    write_string_nodes_csv,
)
from .therapy import (
    CivicTherapyProvider,
    ClinicalTrialsTherapyProvider,
    OncoKBTherapyProvider,
    OpenTargetsTherapyProvider,
    RegulatoryTableTherapyProvider,
    TherapyProvider,
    find_therapy_associations,
)


def _positive_integer(value: str) -> int:
    """Role: Parse and enforce a positive numeric command-line value."""
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def _nonnegative_integer(value: str) -> int:
    """Role: Parse and enforce a non-negative numeric command-line value."""
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be zero or positive")
    return parsed


def _maximum_thirty(value: str) -> int:
    """Role: Compute maximum thirty for the cli workflow."""
    parsed = _positive_integer(value)
    if parsed > 30:
        raise argparse.ArgumentTypeError("must be less than or equal to 30")
    return parsed


def _nonnegative_float(value: str) -> float:
    """Role: Parse and enforce a non-negative numeric command-line value."""
    parsed = float(value)
    if not math.isfinite(parsed) or parsed < 0:
        raise argparse.ArgumentTypeError("must be finite, zero, or positive")
    return parsed


def _evidence_file(value: str) -> tuple[str, Path]:
    """Role: Compute evidence file for the cli workflow."""
    source, separator, path = value.partition("=")
    if not separator or not source.strip() or not path.strip():
        raise argparse.ArgumentTypeError("use SOURCE=PATH, for example clingen=erepo.tsv")
    return source.strip().lower(), Path(path).expanduser()


def build_parser() -> argparse.ArgumentParser:
    """Role: Build the command-line parser and define supported options."""
    parser = argparse.ArgumentParser(
        description=(
            "Find ClinVar variants with ORPHAcodes, detect GDC/TCGA matches, "
            "and retain external evidence in separate layers."
        )
    )
    parser.add_argument("gene", help="HGNC symbol, for example EGFR")
    parser.add_argument(
        "--project",
        action="append",
        default=[],
        help="Specific GDC project; may be repeated (for example, TCGA-LUAD)",
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
        help="ClinVar condition name to exclude; may be repeated",
    )
    parser.add_argument(
        "--classification",
        action="append",
        default=[],
        help=(
            "Allowed ClinVar germline classification; may be repeated. "
            "Use ALL to include VUS, benign, and conflicting classifications."
        ),
    )
    parser.add_argument(
        "--max-clinvar-records",
        type=int,
        help="Limit for smoke tests; omit for a complete analysis",
    )
    parser.add_argument("--email", default=os.getenv("NCBI_EMAIL"))
    parser.add_argument("--api-key", default=os.getenv("NCBI_API_KEY"))
    parser.add_argument(
        "--no-clinvar-details",
        action="store_true",
        help="Skip VCV/RCV/SCV EFetch and use only the ClinVar summary",
    )
    parser.add_argument(
        "--source",
        action="append",
        choices=("alphamissense", "gnomad", "panelapp"),
        default=[],
        help="Additional public source; may be repeated",
    )
    parser.add_argument(
        "--evidence-file",
        action="append",
        type=_evidence_file,
        default=[],
        metavar="SOURCE=PATH",
        help=(
            "CSV/TSV/JSON exported from ClinGen, DECIPHER, LOVD, CIViC, "
            "GenCC, Orphadata, or another source"
        ),
    )
    parser.add_argument(
        "--orphadata-file",
        type=Path,
        help="Official Orphadata en_product6 gene-disease association XML",
    )
    parser.add_argument(
        "--strict-sources",
        action="store_true",
        help="Fail if an additional source cannot be queried or parsed",
    )
    parser.add_argument(
        "--gnomad-max-variants",
        type=_positive_integer,
        default=25,
        help="Maximum candidate alleles queried in gnomAD (maximum 25)",
    )
    parser.add_argument("--output-dir", default="results")
    parser.add_argument(
        "--plots",
        action="store_true",
        help=(
            "Generate lollipop, Venn, recurrence, oncology and rare-disease "
            "heatmaps, and dashboard"
        ),
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
        help="Legacy limit used when specific limits are not provided",
    )
    parser.add_argument(
        "--max-lollipop-residues",
        type=_positive_integer,
        help="Maximum aggregated residues in the lollipop plot (recommended: 80)",
    )
    parser.add_argument(
        "--max-lollipop-labels",
        type=_nonnegative_integer,
        default=15,
        help="Maximum lollipop labels; 0 disables labels",
    )
    parser.add_argument(
        "--max-heatmap-mutations",
        type=_positive_integer,
        help="Independent maximum number of mutations in the heatmap",
    )
    parser.add_argument(
        "--max-rare-disease-mutations",
        type=_maximum_thirty,
        default=30,
        help=(
            "Maximum mutations in the rare-disease heatmap "
            "(1-30; prioritizes cancer-shared mutations)"
        ),
    )
    parser.add_argument(
        "--include-synonymous-plots",
        action="store_true",
        help="Include p.X= changes in the lollipop (omitted by default)",
    )
    parser.add_argument(
        "--string-network",
        action="store_true",
        help=(
            "Retrieve the human STRING v12.0 network, reuse it from SQLite, "
            "and export nodes/edges; --plots also renders it"
        ),
    )
    parser.add_argument(
        "--string-partners",
        type=_nonnegative_integer,
        default=20,
        help="Neighboring proteins added by STRING (0-100; default: 20)",
    )
    parser.add_argument(
        "--string-score",
        type=_nonnegative_integer,
        default=700,
        help="Minimum combined STRING score, 0-1000 (default: 700)",
    )
    parser.add_argument(
        "--string-network-type",
        choices=("functional", "physical"),
        default="functional",
        help="STRING network type (default: functional)",
    )
    parser.add_argument(
        "--structure-inputs",
        action="store_true",
        help=(
            "Generate AlphaFold 3 JSON and Boltz-2 YAML for shared missense "
            "variants"
        ),
    )
    parser.add_argument(
        "--structure-variant",
        action="append",
        default=[],
        metavar="HGVSP",
        help=(
            "Specific protein change, for example L858R; may be repeated and "
            "enables --structure-inputs"
        ),
    )
    parser.add_argument(
        "--structure-max-variants",
        type=_positive_integer,
        default=20,
        help="Maximum generated structural variants (default: 20)",
    )
    parser.add_argument(
        "--protein-sequence-file",
        type=Path,
        help="Local FASTA or canonical sequence; required in offline mode",
    )
    parser.add_argument(
        "--model-seed",
        action="append",
        type=_nonnegative_integer,
        default=[],
        help="AlphaFold 3 seed; may be repeated (default: 1)",
    )
    parser.add_argument(
        "--therapies",
        action="store_true",
        help=(
            "Generate the complementary therapy report with CIViC, "
            "ClinicalTrials.gov, and Open Targets"
        ),
    )
    parser.add_argument(
        "--therapy-source",
        action="append",
        choices=("civic", "oncokb", "clinicaltrials", "opentargets"),
        default=[],
        help=(
            "Therapy source; may be repeated and enables --therapies. "
            "OncoKB requires a license and token"
        ),
    )
    parser.add_argument(
        "--oncokb-token",
        default=os.getenv("ONCOKB_TOKEN"),
        help=(
            "Licensed OncoKB token; define ONCOKB_TOKEN to avoid exposing it "
            "in the shell history"
        ),
    )
    parser.add_argument(
        "--regulatory-file",
        action="append",
        type=_evidence_file,
        default=[],
        metavar="SOURCE=PATH",
        help=(
            "Previously downloaded official FDA/EMA CSV/TSV/JSON; preserves "
            "approval, jurisdiction, and label without scraping"
        ),
    )
    parser.add_argument(
        "--max-therapy-records",
        type=_positive_integer,
        default=100,
        help=(
            "Operational limit per query: ClinicalTrials.gov studies per term, "
            "Open Targets drugs, or CIViC/OncoKB variants"
        ),
    )
    parser.add_argument(
        "--database",
        type=Path,
        help=(
            "Execution SQLite database; defaults to OUTPUT_DIR/rare_oncology.sqlite3"
        ),
    )
    parser.add_argument(
        "--no-database",
        action="store_true",
        help="Disable SQLite history and cache",
    )
    parser.add_argument(
        "--cache-max-age-hours",
        type=_nonnegative_float,
        default=24.0,
        help="Maximum reusable cache age in hours (default: 24)",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Ignore the current cache, query sources, and store another execution",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Use SQLite only; fail if no current compatible execution exists",
    )
    return parser


def _build_therapy_providers(args: argparse.Namespace) -> list[TherapyProvider]:
    """Role: Compute build therapy providers for the cli workflow."""
    if args.therapy_source:
        selected_sources = tuple(dict.fromkeys(args.therapy_source))
    elif args.therapies:
        selected_sources = ("civic", "clinicaltrials", "opentargets")
    else:
        selected_sources = ()

    providers: list[TherapyProvider] = []
    for source in selected_sources:
        if source == "civic":
            providers.append(
                CivicTherapyProvider(max_variants=args.max_therapy_records)
            )
        elif source == "clinicaltrials":
            providers.append(
                ClinicalTrialsTherapyProvider(
                    max_studies=args.max_therapy_records
                )
            )
        elif source == "opentargets":
            providers.append(
                OpenTargetsTherapyProvider(max_drugs=args.max_therapy_records)
            )
        elif source == "oncokb":
            if not args.oncokb_token:
                raise SystemExit(
                    "OncoKB requires ONCOKB_TOKEN or --oncokb-token"
                )
            providers.append(
                OncoKBTherapyProvider(
                    args.oncokb_token,
                    max_variants=args.max_therapy_records,
                )
            )
    for source, path in args.regulatory_file:
        providers.append(
            RegulatoryTableTherapyProvider(
                source,
                path,
                jurisdiction={"fda": "US", "ema": "EU"}.get(source),
            )
        )
    return providers


def _analysis_cache_parameters(
    args: argparse.Namespace,
    allowed_classifications: tuple[str, ...] | None,
) -> dict[str, object]:
    """Role: Compute analysis cache parameters for the cli workflow."""
    return {
        "program": None if args.program.upper() == "ALL" else args.program,
        "projects": sorted(set(args.project)),
        "min_cancer_patients": args.min_patients,
        "max_clinvar_records": args.max_clinvar_records,
        "exclude_conditions": sorted(set(args.exclude_condition)),
        "allowed_germline_classifications": (
            None
            if allowed_classifications is None
            else sorted(set(allowed_classifications))
        ),
        "include_clinvar_details": not args.no_clinvar_details,
        "sources": list(dict.fromkeys(args.source)),
        "evidence_files": [
            {"source": source, "path": path}
            for source, path in args.evidence_file
        ],
        "orphadata_file": args.orphadata_file,
        "gnomad_max_variants": args.gnomad_max_variants,
        "strict_sources": args.strict_sources,
    }


def _therapy_cache_parameters(args: argparse.Namespace) -> dict[str, object]:
    """Role: Compute therapy cache parameters for the cli workflow."""
    if args.therapy_source:
        sources = list(dict.fromkeys(args.therapy_source))
    elif args.therapies:
        sources = ["civic", "clinicaltrials", "opentargets"]
    else:
        sources = []
    return {
        "sources": sources,
        "regulatory_files": [
            {"source": source, "path": path}
            for source, path in args.regulatory_file
        ],
        "max_therapy_records": args.max_therapy_records,
        "strict_sources": args.strict_sources,
    }


def main(argv: list[str] | None = None) -> int:
    """Role: Run the module entry point and coordinate its user-facing workflow."""
    args = build_parser().parse_args(argv)
    if args.gnomad_max_variants > 25:
        raise SystemExit("--gnomad-max-variants cannot exceed 25")
    if args.string_partners > 100:
        raise SystemExit("--string-partners cannot exceed 100")
    if args.string_score > 1000:
        raise SystemExit("--string-score cannot exceed 1000")
    if args.offline and args.refresh:
        raise SystemExit("--offline and --refresh cannot be used together")
    if args.no_database and (args.offline or args.refresh or args.database):
        raise SystemExit(
            "--no-database cannot be combined with --database, --offline, or --refresh"
        )
    output_dir = Path(args.output_dir)
    database_path = args.database or output_dir / "rare_oncology.sqlite3"
    clinvar = NcbiClinVarClient(
        email=args.email,
        api_key=args.api_key,
        include_details=not args.no_clinvar_details,
    )
    providers = []
    for source in dict.fromkeys(args.source):
        if source == "alphamissense":
            providers.append(
                AlphaMissenseEvidenceProvider(
                    database=None if args.no_database else database_path
                )
            )
        elif source == "gnomad":
            providers.append(
                GnomadEvidenceProvider(max_variants=args.gnomad_max_variants)
            )
        elif source == "panelapp":
            providers.append(PanelAppEvidenceProvider())
    providers.extend(
        LocalTableEvidenceProvider(source, path)
        for source, path in args.evidence_file
    )
    if args.orphadata_file:
        providers.append(OrphadataXmlEvidenceProvider(args.orphadata_file))
    allowed_classifications = (
        None
        if any(value.upper() == "ALL" for value in args.classification)
        else tuple(args.classification) or DEFAULT_ALLOWED_CLASSIFICATIONS
    )
    analysis_kwargs = {
        "min_cancer_patients": args.min_patients,
        "program": None if args.program.upper() == "ALL" else args.program,
        "projects": args.project,
        "max_clinvar_records": args.max_clinvar_records,
        "exclude_conditions": args.exclude_condition,
        "allowed_germline_classifications": allowed_classifications,
        "clinvar_client": clinvar,
        "evidence_providers": providers,
        "strict_evidence_sources": args.strict_sources,
    }
    analysis_cache = None
    try:
        if args.no_database:
            result = find_shared_mutations(args.gene, **analysis_kwargs)
        else:
            analysis_cache = find_shared_mutations_cached(
                args.gene,
                database=database_path,
                cache_parameters=_analysis_cache_parameters(
                    args, allowed_classifications
                ),
                max_age_hours=args.cache_max_age_hours,
                refresh=args.refresh,
                offline=args.offline,
                **analysis_kwargs,
            )
            result = analysis_cache.result
    except CacheMissError as exc:
        raise SystemExit(str(exc)) from exc

    stem = result.gene.lower()
    json_path = write_json(result, output_dir / f"{stem}_report.json")
    csv_path = write_csv(result, output_dir / f"{stem}_shared.csv")
    evidence_path = write_evidence_csv(
        result, output_dir / f"{stem}_evidence.csv"
    )
    print(
        f"{result.gene}: {result.rare_variant_count} rare variants, "
        f"{result.cancer_observation_count} oncology observations, "
        f"{result.shared_mutation_count} candidate matches, "
        f"{result.evidence_record_count} evidence records."
    )
    print(f"JSON: {json_path}")
    print(f"CSV:  {csv_path}")
    print(f"Evidence CSV: {evidence_path}")
    if analysis_cache is not None:
        print(
            f"SQLite molecular: {'HIT' if analysis_cache.cache_hit else 'WEB→DB'} "
            f"execution {analysis_cache.execution_id} "
            f"[{analysis_cache.database_path}]"
        )
    for status in result.source_statuses:
        suffix = f" ({status.message})" if status.message else ""
        print(
            f"Source {status.source}: {status.status}, "
            f"{status.record_count} records [{status.evidence_layer}]{suffix}"
        )
    string_network_run = None
    if args.string_network:
        try:
            string_network_run = get_string_network(
                result.gene,
                database=None if args.no_database else database_path,
                required_score=args.string_score,
                network_type=args.string_network_type,
                added_nodes=args.string_partners,
                refresh=args.refresh,
                offline=args.offline,
            )
        except (StringNetworkError, ValueError) as exc:
            raise SystemExit(str(exc)) from exc
        network = string_network_run.network
        string_json = write_string_network_json(
            network, output_dir / f"{stem}_string_network.json"
        )
        string_nodes = write_string_nodes_csv(
            network, output_dir / f"{stem}_string_nodes.csv"
        )
        string_edges = write_string_edges_csv(
            network, output_dir / f"{stem}_string_edges.csv"
        )
        string_status = (
            "HIT SQLite"
            if string_network_run.cache_hit
            else "WEB→SQLite"
            if string_network_run.database_path
            else "WEB"
        )
        print(
            f"STRING v{network.string_version}: {len(network.nodes)} nodes, "
            f"{len(network.edges)} edges [{string_status}]."
        )
        print(f"STRING JSON:   {string_json}")
        print(f"STRING nodes:  {string_nodes}")
        print(f"STRING edges: {string_edges}")
    therapy_requested = bool(
        args.therapies or args.therapy_source or args.regulatory_file
    )
    if therapy_requested:
        therapy_providers = _build_therapy_providers(args)
        therapy_cache = None
        try:
            if args.no_database:
                therapy_result = find_therapy_associations(
                    result,
                    providers=therapy_providers,
                    strict=args.strict_sources,
                )
            else:
                therapy_cache = find_therapy_associations_cached(
                    result,
                    database=database_path,
                    providers=therapy_providers,
                    cache_parameters=_therapy_cache_parameters(args),
                    max_age_hours=args.cache_max_age_hours,
                    refresh=args.refresh,
                    offline=args.offline,
                    strict=args.strict_sources,
                )
                therapy_result = therapy_cache.result
        except CacheMissError as exc:
            raise SystemExit(str(exc)) from exc
        therapy_json = write_therapy_json(
            therapy_result, output_dir / f"{stem}_therapies.json"
        )
        therapy_csv = write_therapy_csv(
            therapy_result, output_dir / f"{stem}_therapies.csv"
        )
        print(
            f"Therapies: {therapy_result.association_count} associations "
            "across three independent layers."
        )
        print(f"Therapies JSON: {therapy_json}")
        print(f"Therapies CSV:  {therapy_csv}")
        if therapy_cache is not None:
            print(
                f"Therapy SQLite: {'HIT' if therapy_cache.cache_hit else 'WEB→DB'} "
                f"execution {therapy_cache.execution_id} "
                f"[{therapy_cache.database_path}]"
            )
        for status in therapy_result.source_statuses:
            suffix = f" ({status.message})" if status.message else ""
            print(
                f"Therapy source {status.source}: {status.status}, "
                f"{status.record_count} records [{status.evidence_layer}]"
                f"{suffix}"
            )
    structure_requested = bool(
        args.structure_inputs or args.structure_variant or args.protein_sequence_file
    )
    if structure_requested:
        try:
            sequence = (
                read_protein_sequence(args.protein_sequence_file)
                if args.protein_sequence_file is not None
                else None
            )
            structure_result = generate_variant_structure_inputs(
                result,
                output_dir / "structures",
                variants=tuple(args.structure_variant) or None,
                protein_sequence=sequence,
                sequence_database=None if args.no_database else database_path,
                offline=args.offline,
                model_seeds=tuple(args.model_seed) or (1,),
                max_variants=args.structure_max_variants,
            )
        except StructureInputError as exc:
            raise SystemExit(str(exc)) from exc
        print(
            f"Structures: {len(structure_result.artifacts)} files, "
            f"{len(structure_result.skipped)} variants skipped."
        )
        print(f"Manifest: {structure_result.manifest_path}")
        for artifact in structure_result.artifacts:
            print(
                f"{artifact.platform} p.{artifact.protein_change}: "
                f"{artifact.path}"
            )
        for skipped in structure_result.skipped:
            print(f"Skipped p.{skipped.protein_change}: {skipped.reason}")
    if args.plots:
        offline_annotation = None
        if args.offline:
            positions = [
                record.amino_acid_position
                for record in result.mutation_plot_data
                if record.amino_acid_position is not None
            ]
            offline_annotation = inferred_annotation(result.gene, positions)
        try:
            plot_paths = generate_gene_plots(
                result,
                output_dir / "plots",
                formats=tuple(args.plot_format) or ("png", "svg"),
                max_lollipop_mutations=(
                    args.max_lollipop_residues or args.max_plot_mutations
                ),
                max_lollipop_labels=args.max_lollipop_labels,
                max_heatmap_mutations=(
                    args.max_heatmap_mutations or args.max_plot_mutations
                ),
                max_rare_disease_mutations=args.max_rare_disease_mutations,
                include_synonymous_lollipop=args.include_synonymous_plots,
                protein_annotation=offline_annotation,
            )
        except PlotDependencyError as exc:
            raise SystemExit(str(exc)) from exc
        for name, paths in plot_paths.items():
            print(f"{name.capitalize()}: {', '.join(str(path) for path in paths)}")
        if string_network_run is not None:
            try:
                string_paths = generate_string_network_plots(
                    string_network_run.network,
                    output_dir / "plots",
                    formats=tuple(args.plot_format) or ("png", "svg"),
                )
            except PlotDependencyError as exc:
                raise SystemExit(str(exc)) from exc
            print(
                "STRING network: "
                + ", ".join(str(path) for path in string_paths)
            )
    return 0
