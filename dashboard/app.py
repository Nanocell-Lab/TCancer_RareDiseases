"""Streamlit entry point for the rare-oncology 2.1 dashboard."""

from __future__ import annotations

import hashlib
import os
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Mapping


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
SRC = HERE / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import streamlit as st  # noqa: E402

from rare_oncology_dashboard.catalog import (  # noqa: E402
    ANALYSIS_KIND,
    THERAPY_KIND,
    ExecutionSnapshot,
    SQLiteCatalog,
    StringNetworkSnapshot,
    discover_databases,
)
from rare_oncology_dashboard.export import (  # noqa: E402
    files_to_zip,
    payload_to_json,
    records_to_csv,
)
from rare_oncology_dashboard.figures import (  # noqa: E402
    alphamissense_figure,
    heatmap_figure,
    lollipop_figure,
    overlap_figure,
    recurrence_figure,
    residue_displacement_figure,
    string_evidence_heatmap_figure,
    string_network_figure,
    string_partner_scores_figure,
    structure_overlay_figure,
    therapy_summary_figure,
)
from rare_oncology_dashboard.runner import (  # noqa: E402
    AnalysisRunner,
    RunOptions,
    RunnerError,
    StringRunOptions,
    StructureRunOptions,
    normalize_gene,
)
from rare_oncology_dashboard.structure_compare import (  # noqa: E402
    PDBComparisonError,
    compare_pdb_models,
)
from rare_oncology_dashboard.transforms import (  # noqa: E402
    alphamissense_rows,
    lollipop_rows,
    oncology_heatmap,
    overlap_counts,
    rare_disease_heatmap,
    shared_mutation_rows,
    source_status_rows,
    string_graph_data,
    string_partner_rows,
    structure_variant_rows,
    therapy_rows,
)


RESULTS_ROOT = PROJECT_ROOT / "2.1" / "results"
DEFAULT_MASTER = RESULTS_ROOT / "dashboard_master.sqlite3"


def main() -> None:
    """Role: Run the module entry point and coordinate its user-facing workflow."""
    st.set_page_config(
        page_title="Rare × Oncology Gene Explorer",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    _style()

    master_text = st.sidebar.text_input(
        "Writable master SQLite",
        value=str(DEFAULT_MASTER),
        help="New 2.1 executions are written here. Existing databases remain read-only.",
    )
    master_path = Path(master_text).expanduser().resolve()
    discovered = discover_databases(RESULTS_ROOT, master_path)
    selected_paths = _database_selector(discovered, master_path)
    catalog = SQLiteCatalog(tuple(selected_paths))
    entries = catalog.genes()
    entry_by_gene = {entry.gene: entry for entry in entries}
    genes = list(entry_by_gene)

    st.sidebar.markdown("---")
    st.sidebar.subheader("Gene selection")
    preferred_index = genes.index("TGFBR1") if "TGFBR1" in genes else 0
    selected_gene = st.sidebar.selectbox(
        "Stored genes",
        options=genes or [""],
        index=preferred_index if genes else 0,
        format_func=lambda gene: _gene_label(gene, entry_by_gene),
    )
    custom_gene = st.sidebar.text_input(
        "Or enter an HGNC symbol",
        placeholder="e.g. EGFR, TGFBR1, ACVR1",
    )
    requested_gene = custom_gene.strip() or selected_gene
    try:
        gene = normalize_gene(requested_gene) if requested_gene else ""
    except ValueError as exc:
        st.sidebar.error(str(exc))
        gene = ""

    st.sidebar.markdown("---")
    with st.sidebar.expander("Run options for missing or refreshed genes"):
        evidence_sources = tuple(
            st.multiselect(
                "Additional evidence",
                ["alphamissense", "gnomad", "panelapp"],
                default=["alphamissense", "gnomad", "panelapp"],
            )
        )
        include_therapies = st.checkbox("Include therapy layers", value=True)
        oncokb_available = bool(os.environ.get("ONCOKB_TOKEN"))
        include_oncokb = st.checkbox(
            "Include licensed OncoKB",
            value=False,
            disabled=not oncokb_available,
            help=(
                "Uses ONCOKB_TOKEN from the environment; the token is never placed "
                "in the command line or dashboard."
            ),
        )
        max_clinvar = int(
            st.number_input(
                "ClinVar record limit (0 = complete)",
                min_value=0,
                max_value=100000,
                value=0,
                step=100,
            )
        )
        max_therapy = int(
            st.number_input(
                "Maximum therapy records per query",
                min_value=1,
                max_value=1000,
                value=100,
            )
        )
        include_string = st.checkbox("Include STRING partner network", value=True)
        string_score = int(
            st.number_input(
                "Minimum STRING score",
                min_value=0,
                max_value=1000,
                value=700,
                step=50,
            )
        )
        string_partners = int(
            st.number_input(
                "STRING partner limit",
                min_value=0,
                max_value=100,
                value=20,
                step=5,
            )
        )
        string_network_type = st.selectbox(
            "STRING network type", ("functional", "physical")
        )

    analysis_snapshot = catalog.latest(gene, ANALYSIS_KIND) if gene else None
    string_snapshot = catalog.latest_string_network(gene) if gene else None
    col_run, col_refresh = st.sidebar.columns(2)
    run_missing = col_run.button(
        "Run missing",
        type="primary",
        width="stretch",
        disabled=not gene or analysis_snapshot is not None,
    )
    refresh = col_refresh.button(
        "Refresh web",
        width="stretch",
        disabled=not gene,
    )
    if run_missing or refresh:
        _run_analysis(
            gene,
            master_path,
            RunOptions(
                include_therapies=include_therapies,
                include_oncokb=include_oncokb,
                evidence_sources=evidence_sources,
                max_clinvar_records=max_clinvar or None,
                max_therapy_records=max_therapy,
                include_string_network=include_string,
                string_score=string_score,
                string_partners=string_partners,
                string_network_type=string_network_type,
                force_refresh=refresh,
            ),
        )
    string_action = st.sidebar.button(
        "Fetch STRING" if string_snapshot is None else "Refresh STRING",
        width="stretch",
        disabled=not gene,
        help="Runs only the cache-first STRING connector; it does not repeat ClinVar or GDC.",
    )
    if string_action:
        _run_string_network(
            gene,
            master_path,
            StringRunOptions(
                required_score=string_score,
                added_nodes=string_partners,
                network_type=string_network_type,
                force_refresh=string_snapshot is not None,
            ),
        )

    _header(len(genes), len(selected_paths), master_path)
    for warning in catalog.warnings():
        st.warning(warning)
    if last_run := st.session_state.pop("last_run", None):
        with st.expander("Most recent 2.1 execution", expanded=False):
            st.success(f"Completed {last_run['gene']} into {last_run['database']}")
            st.code(last_run["command"], language="bash")
            if last_run["stdout"]:
                st.text(last_run["stdout"])

    if not gene:
        st.info("Choose a stored gene or enter an HGNC symbol to begin.")
        return
    if analysis_snapshot is None:
        st.warning(
            f"{gene} is not present in the selected SQLite databases. "
            "Press **Run missing** to execute the cache-first 2.1 pipeline and populate the master database."
        )
        _catalog_overview(entries)
        return

    analysis_snapshot = _snapshot_selector(catalog, gene, ANALYSIS_KIND, analysis_snapshot)
    therapy_snapshot = catalog.latest(gene, THERAPY_KIND)
    if therapy_snapshot is not None:
        therapy_snapshot = _snapshot_selector(
            catalog, gene, THERAPY_KIND, therapy_snapshot, sidebar=False
        )
    if string_snapshot is not None:
        string_snapshot = _string_snapshot_selector(catalog, gene, string_snapshot)
    _render_gene(analysis_snapshot, therapy_snapshot, string_snapshot, master_path)


def _database_selector(
    discovered: tuple[Path, ...], master_path: Path
) -> tuple[Path, ...]:
    """Role: Compute database selector for the app workflow."""
    st.sidebar.subheader("2.1 data sources")
    labels = {str(path): path for path in discovered}
    selected = st.sidebar.multiselect(
        "SQLite databases",
        options=list(labels),
        default=list(labels),
        format_func=lambda value: _relative(Path(value)),
        help="The dashboard federates these databases without copying their payloads.",
    )
    selected_paths = {labels[value] for value in selected}
    if master_path.is_file():
        selected_paths.add(master_path.resolve())
    st.sidebar.caption(f"{len(selected_paths)} database(s) selected (master always linked)")
    return tuple(sorted(selected_paths, key=str))


def _snapshot_selector(
    catalog: SQLiteCatalog,
    gene: str,
    kind: str,
    latest: ExecutionSnapshot,
    *,
    sidebar: bool = True,
) -> ExecutionSnapshot:
    """Role: Compute snapshot selector for the app workflow."""
    snapshots = catalog.snapshots(gene, kind)
    if len(snapshots) == 1:
        return latest
    location = st.sidebar if sidebar else st
    label = "Molecular snapshot" if kind == ANALYSIS_KIND else "Therapy snapshot"
    selected_identity = location.selectbox(
        label,
        options=[snapshot.identity for snapshot in snapshots],
        key=f"snapshot-{kind}-{gene}",
    )
    return next(snapshot for snapshot in snapshots if snapshot.identity == selected_identity)


def _string_snapshot_selector(
    catalog: SQLiteCatalog,
    gene: str,
    latest: StringNetworkSnapshot,
) -> StringNetworkSnapshot:
    """Role: Compute string snapshot selector for the app workflow."""
    snapshots = catalog.string_networks(gene)
    if len(snapshots) == 1:
        return latest
    selected_identity = st.selectbox(
        "STRING snapshot",
        options=[snapshot.identity for snapshot in snapshots],
        key=f"snapshot-string-{gene}",
    )
    return next(snapshot for snapshot in snapshots if snapshot.identity == selected_identity)


def _run_analysis(gene: str, master_path: Path, options: RunOptions) -> None:
    """Role: Compute run analysis for the app workflow."""
    runner = AnalysisRunner(PROJECT_ROOT, master_path)
    with st.spinner(
        f"Running rare-oncology 2.1 for {gene}. SQLite is checked by 2.1 before web sources..."
    ):
        try:
            outcome = runner.run(gene, options)
        except (RunnerError, ValueError) as exc:
            st.error(str(exc))
            return
    st.session_state["last_run"] = {
        "gene": outcome.gene,
        "database": str(outcome.master_database),
        "command": outcome.command_display,
        "stdout": outcome.stdout,
    }
    st.rerun()


def _run_string_network(
    gene: str, master_path: Path, options: StringRunOptions
) -> None:
    """Role: Compute run string network for the app workflow."""
    runner = AnalysisRunner(PROJECT_ROOT, master_path)
    with st.spinner(
        f"Checking SQLite and retrieving missing STRING partners for {gene}..."
    ):
        try:
            outcome = runner.run_string_network(gene, options)
        except (RunnerError, ValueError) as exc:
            st.error(str(exc))
            return
    st.session_state["last_run"] = {
        "gene": outcome.gene,
        "database": str(outcome.master_database),
        "command": outcome.command_display,
        "stdout": outcome.stdout,
    }
    st.rerun()


def _render_gene(
    analysis_snapshot: ExecutionSnapshot,
    therapy_snapshot: ExecutionSnapshot | None,
    string_snapshot: StringNetworkSnapshot | None,
    master_path: Path,
) -> None:
    """Role: Render the gene section in the interactive dashboard."""
    result = analysis_snapshot.result
    therapy_result = therapy_snapshot.result if therapy_snapshot else None
    gene = analysis_snapshot.gene
    st.markdown(f"## {gene} molecular landscape")
    st.caption(
        f"Stored execution #{analysis_snapshot.execution_id} · "
        f"{analysis_snapshot.created_at} · {_relative(analysis_snapshot.database_path)}"
    )
    string_rows = string_partner_rows(string_snapshot.result if string_snapshot else None)
    cols = st.columns(6)
    cols[0].metric("Rare-disease variants", _count(result, "rare_variant_count"))
    cols[1].metric("Cancer observations", _count(result, "cancer_observation_count"))
    cols[2].metric("Shared mutations", _count(result, "shared_mutation_count"))
    cols[3].metric("Evidence records", _count(result, "evidence_record_count"))
    cols[4].metric(
        "Therapy associations",
        _count(therapy_result, "association_count") if therapy_result else "—",
    )
    cols[5].metric("STRING direct partners", len(string_rows) if string_snapshot else "—")
    st.info(str(result.get("caveat", "Research-use result; confirm variant identity clinically.")))

    (
        overview_tab,
        mutation_tab,
        rare_tab,
        alpha_tab,
        therapy_tab,
        string_tab,
        structure_tab,
        provenance_tab,
    ) = st.tabs(
        [
            "Overview",
            "Mutations",
            "Rare diseases",
            "AlphaMissense",
            "Therapies",
            "STRING partners",
            "Structures",
            "Provenance",
        ]
    )
    shared_rows = shared_mutation_rows(result)
    with overview_tab:
        left, right = st.columns([1, 2])
        with left:
            st.plotly_chart(overlap_figure(overlap_counts(result)), width="stretch")
        with right:
            st.plotly_chart(recurrence_figure(shared_rows, top=20), width="stretch")
        _download_table("Download shared mutations", shared_rows, f"{gene.lower()}_shared.csv")
        st.dataframe(shared_rows, width="stretch", hide_index=True)

    with mutation_tab:
        max_residues = st.slider(
            "Maximum lollipop residues",
            min_value=20,
            max_value=200,
            value=100,
            step=10,
        )
        lollipops = lollipop_rows(result, max_residues=max_residues)
        st.plotly_chart(lollipop_figure(lollipops, gene), width="stretch")
        st.caption(
            "Stem height uses log2(cancer patients + 1). Diamonds/bold borders indicate "
            "genomic-and-protein confirmation; hover exposes diseases and projects."
        )
        cancer_matrix = oncology_heatmap(result, top=30)
        st.plotly_chart(
            heatmap_figure(
                cancer_matrix,
                row_key="mutations",
                column_key="projects",
                title="Top oncology mutations by project (★ shared with rare disease)",
                colorbar_title="Patients",
            ),
            width="stretch",
        )

    with rare_tab:
        rare_matrix = rare_disease_heatmap(result, top=30)
        st.plotly_chart(
            heatmap_figure(
                rare_matrix,
                row_key="diseases",
                column_key="mutations",
                title="Top 30 rare-disease mutation associations (★ shared with cancer)",
                colorbar_title="ClinVar associations",
            ),
            width="stretch",
        )
        st.caption(
            "Shared cancer mutations are selected before rare-only mutations; within each "
            "group, disease-association abundance determines rank."
        )

    with alpha_tab:
        alpha_rows = alphamissense_rows(result)
        st.plotly_chart(alphamissense_figure(alpha_rows, gene), width="stretch")
        _download_table(
            "Download AlphaMissense table", alpha_rows, f"{gene.lower()}_alphamissense.csv"
        )
        st.dataframe(alpha_rows, width="stretch", hide_index=True)

    with therapy_tab:
        therapies = therapy_rows(therapy_result)
        if therapy_snapshot:
            st.caption(
                f"Therapy execution #{therapy_snapshot.execution_id} · "
                f"{therapy_snapshot.created_at} · {_relative(therapy_snapshot.database_path)}"
            )
        st.plotly_chart(therapy_summary_figure(therapies), width="stretch")
        if therapy_result:
            st.warning(str(therapy_result.get("caveat", "Therapy evidence is not a clinical recommendation.")))
        layer_options = sorted({str(row["layer"]) for row in therapies})
        selected_layers = st.multiselect(
            "Therapy layers",
            layer_options,
            default=layer_options,
        )
        filtered = [row for row in therapies if row["layer"] in selected_layers]
        _download_table("Download therapy table", filtered, f"{gene.lower()}_therapies.csv")
        st.dataframe(filtered, width="stretch", hide_index=True)

    with string_tab:
        if string_snapshot is None:
            st.warning(
                "No STRING snapshot is stored for this gene. Use **Fetch STRING** "
                "in the sidebar; ClinVar and GDC will not be rerun."
            )
        else:
            string_result = string_snapshot.result
            st.caption(
                f"STRING v{string_snapshot.string_version} · Homo sapiens "
                f"({string_snapshot.species}) · {string_snapshot.network_type} · "
                f"minimum score {string_snapshot.required_score} · "
                f"{string_snapshot.created_at} · {_relative(string_snapshot.database_path)}"
            )
            graph = string_graph_data(string_result)
            st.plotly_chart(string_network_figure(graph, gene), width="stretch")
            st.caption(
                f"Showing {graph['displayed_edges']} of {graph['total_edges']} stored edges. "
                "All direct query–partner associations are retained; only secondary edges are capped."
            )
            score_col, evidence_col = st.columns(2)
            with score_col:
                st.plotly_chart(
                    string_partner_scores_figure(string_rows), width="stretch"
                )
            with evidence_col:
                st.plotly_chart(
                    string_evidence_heatmap_figure(string_rows), width="stretch"
                )
            _download_table(
                "Download STRING partners",
                string_rows,
                f"{gene.lower()}_string_partners.csv",
            )
            st.dataframe(string_rows, width="stretch", hide_index=True)
            st.warning(
                str(
                    string_result.get(
                        "caveat",
                        "STRING associations do not by themselves prove physical binding or causality.",
                    )
                )
            )
            if source_url := string_result.get("source_url"):
                st.markdown(f"[Open this network in STRING]({source_url})")

    with structure_tab:
        _render_structures(
            gene,
            result,
            analysis_snapshot.execution_id,
            master_path,
        )

    with provenance_tab:
        statuses = source_status_rows(result, therapy_result)
        st.subheader("Source status")
        st.dataframe(statuses, width="stretch", hide_index=True)
        left, right = st.columns(2)
        with left:
            st.subheader("Molecular request")
            st.json(analysis_snapshot.parameters)
        with right:
            st.subheader("Stored source URLs")
            st.write(list(result.get("sources", [])))
        download_count = 1 + int(bool(therapy_result)) + int(bool(string_snapshot))
        downloads = st.columns(download_count)
        downloads[0].download_button(
            "Download molecular JSON",
            data=payload_to_json(result),
            file_name=f"{gene.lower()}_report.json",
            mime="application/json",
        )
        if therapy_result:
            downloads[1].download_button(
                "Download therapy JSON",
                data=payload_to_json(therapy_result),
                file_name=f"{gene.lower()}_therapies.json",
                mime="application/json",
            )
        if string_snapshot:
            downloads[-1].download_button(
                "Download STRING JSON",
                data=payload_to_json(string_snapshot.result),
                file_name=f"{gene.lower()}_string_network.json",
                mime="application/json",
            )


def _render_structures(
    gene: str,
    result: Mapping[str, Any],
    execution_id: int,
    master_path: Path,
) -> None:
    """Role: Render the structures section in the interactive dashboard."""
    input_tab, comparison_tab = st.tabs(
        ["AlphaFold 3 / Boltz-2 inputs", "PDB wild-type vs variant"]
    )
    with input_tab:
        candidates = structure_variant_rows(result)
        by_change = {str(row["protein_change"]): row for row in candidates}
        options = list(by_change)
        defaults = [
            change
            for change in options
            if by_change[change]["category"] == "Shared"
        ][:3]
        selected = st.multiselect(
            "Missense variants",
            options,
            default=defaults,
            format_func=lambda change: (
                f"p.{change} · {by_change[change]['category']} · "
                f"{by_change[change]['cancer_patient_count']} cancer patient(s)"
            ),
            key=f"structure-variants-{gene}-{execution_id}",
        )
        sequence_upload = st.file_uploader(
            "Optional canonical protein FASTA",
            type=("fa", "fasta", "faa", "txt"),
            help=(
                "If omitted, 2.1 checks protein_sequence_cache first and only then "
                "queries reviewed UniProt. The reference residue must match."
            ),
            key=f"structure-sequence-{gene}-{execution_id}",
        )
        seed_text = st.text_input(
            "AlphaFold 3 model seeds",
            value="1",
            help="Comma-separated non-negative integers, for example 1,2,3.",
            key=f"structure-seeds-{gene}-{execution_id}",
        )
        offline_sequence = st.checkbox(
            "SQLite/FASTA only (do not query UniProt)",
            value=False,
            key=f"structure-offline-{gene}-{execution_id}",
        )
        generate = st.button(
            "Generate prediction inputs",
            type="primary",
            disabled=not selected,
            key=f"generate-structures-{gene}-{execution_id}",
        )
        bundle_key = f"structure-bundle-{gene}-{execution_id}"
        sequence_bytes = (
            sequence_upload.getvalue() if sequence_upload is not None else b""
        )
        request_signature = (
            tuple(selected),
            seed_text,
            offline_sequence,
            hashlib.sha256(sequence_bytes).hexdigest(),
        )
        if generate:
            try:
                seeds = _parse_model_seeds(seed_text)
                sequence = (
                    sequence_bytes.decode("utf-8")
                    if sequence_upload is not None
                    else None
                )
                runner = AnalysisRunner(PROJECT_ROOT, master_path)
                with st.spinner(
                    "Validating the canonical sequence and building mutant inputs..."
                ):
                    outcome = runner.run_structure_inputs(
                        gene,
                        result,
                        StructureRunOptions(
                            variants=tuple(selected),
                            model_seeds=seeds,
                            protein_sequence=sequence,
                            offline=offline_sequence,
                        ),
                    )
                st.session_state[bundle_key] = {
                    "request_signature": request_signature,
                    "outcome": outcome,
                }
            except (RunnerError, UnicodeDecodeError, ValueError) as exc:
                st.error(str(exc))
        stored_bundle = st.session_state.get(bundle_key)
        outcome = (
            stored_bundle["outcome"]
            if stored_bundle is not None
            and stored_bundle["request_signature"] == request_signature
            else None
        )
        if outcome is not None:
            st.success(
                f"Generated {len(outcome.files) - 1} platform inputs for {gene}."
            )
            st.caption(
                "Generated by the focused 2.1 structure workflow after canonical "
                "reference-residue validation."
            )
            archive = files_to_zip(
                [(item.filename, item.content) for item in outcome.files]
            )
            st.download_button(
                "Download all structure inputs (.zip)",
                data=archive,
                file_name=f"{gene.lower()}_structure_inputs.zip",
                mime="application/zip",
                key=f"download-structure-zip-{gene}-{execution_id}",
            )
            for item in outcome.files:
                label = (
                    f"Download {item.platform} p.{item.protein_change}"
                    if item.protein_change
                    else "Download generation manifest"
                )
                st.download_button(
                    label,
                    data=item.content,
                    file_name=item.filename,
                    mime=item.mime_type,
                    key=f"download-{item.filename}-{execution_id}",
                )
            if outcome.skipped:
                st.warning("Some selected variants were not generated.")
                st.dataframe(list(outcome.skipped), width="stretch", hide_index=True)
        if not candidates:
            st.info(
                "This molecular snapshot has no simple missense variants eligible "
                "for sequence substitution."
            )
        st.caption(
            "AlphaFold 3 receives JSON dialect version 4. Boltz-2 receives its "
            "native YAML schema. Inputs model one canonical protein chain and do "
            "not automatically add ligands, glycans, membrane or oligomeric partners."
        )

    with comparison_tab:
        st.write(
            "Upload both models. The dashboard matches Cα atoms by chain and residue "
            "identifier, globally aligns the variant with Kabsch, and reports shifts."
        )
        wildtype_upload = st.file_uploader(
            "Wild-type PDB",
            type=("pdb",),
            key=f"wildtype-pdb-{gene}-{execution_id}",
        )
        variant_upload = st.file_uploader(
            "Variant PDB",
            type=("pdb",),
            key=f"variant-pdb-{gene}-{execution_id}",
        )
        compare = st.button(
            "Compare PDB models",
            disabled=wildtype_upload is None or variant_upload is None,
            key=f"compare-pdb-{gene}-{execution_id}",
        )
        comparison_key = f"pdb-comparison-{gene}-{execution_id}"
        if compare and wildtype_upload is not None and variant_upload is not None:
            wildtype_data = wildtype_upload.getvalue()
            variant_data = variant_upload.getvalue()
            try:
                comparison = compare_pdb_models(
                    wildtype_data,
                    variant_data,
                    wildtype_name=wildtype_upload.name,
                    variant_name=variant_upload.name,
                )
            except PDBComparisonError as exc:
                st.error(str(exc))
            else:
                st.session_state[comparison_key] = {
                    "wildtype_sha256": hashlib.sha256(wildtype_data).hexdigest(),
                    "variant_sha256": hashlib.sha256(variant_data).hexdigest(),
                    "result": comparison,
                }
        stored = st.session_state.get(comparison_key)
        if stored is not None and wildtype_upload is not None and variant_upload is not None:
            current_hashes = (
                hashlib.sha256(wildtype_upload.getvalue()).hexdigest(),
                hashlib.sha256(variant_upload.getvalue()).hexdigest(),
            )
            if current_hashes == (
                stored["wildtype_sha256"],
                stored["variant_sha256"],
            ):
                _render_pdb_comparison(gene, stored["result"], execution_id)


def _render_pdb_comparison(
    gene: str, comparison: Mapping[str, Any], execution_id: int
) -> None:
    """Role: Render the pdb comparison section in the interactive dashboard."""
    metrics = st.columns(4)
    metrics[0].metric(
        "Global Cα RMSD",
        f"{float(comparison['global_ca_rmsd_angstrom']):.3f} Å",
    )
    metrics[1].metric("Common Cα atoms", int(comparison["common_ca_count"]))
    metrics[2].metric(
        "Coverage", f"{100 * float(comparison['coverage_fraction']):.1f}%"
    )
    metrics[3].metric(
        "Maximum displacement",
        f"{float(comparison['max_ca_displacement_angstrom']):.3f} Å",
    )
    st.plotly_chart(structure_overlay_figure(comparison, gene), width="stretch")
    st.plotly_chart(residue_displacement_figure(comparison), width="stretch")
    residues = list(comparison.get("residues", []))
    st.dataframe(residues, width="stretch", hide_index=True)
    downloads = st.columns(2)
    downloads[0].download_button(
        "Download comparison JSON",
        data=payload_to_json(comparison),
        file_name=f"{gene.lower()}_pdb_comparison.json",
        mime="application/json",
        key=f"download-pdb-json-{gene}-{execution_id}",
    )
    downloads[1].download_button(
        "Download residue shifts CSV",
        data=records_to_csv(residues),
        file_name=f"{gene.lower()}_pdb_residue_shifts.csv",
        mime="text/csv",
        key=f"download-pdb-csv-{gene}-{execution_id}",
    )
    st.warning(str(comparison.get("caveat", "Interpret structural differences cautiously.")))


def _parse_model_seeds(value: str) -> tuple[int, ...]:
    """Role: Parse and validate model seeds from its source representation."""
    try:
        seeds = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as exc:
        raise ValueError("Model seeds must be comma-separated integers") from exc
    if not seeds or any(seed < 0 for seed in seeds):
        raise ValueError("Model seeds must contain non-negative integers")
    return tuple(dict.fromkeys(seeds))


def _catalog_overview(entries: tuple[Any, ...]) -> None:
    """Role: Compute catalog overview for the app workflow."""
    st.subheader("Stored gene catalog")
    st.dataframe(
        [
            {
                "gene": entry.gene,
                "molecular_runs": entry.analysis_runs,
                "therapy_runs": entry.therapy_runs,
                "latest_molecular": entry.latest_analysis_at,
                "latest_therapy": entry.latest_therapy_at,
                "databases": ", ".join(_relative(path) for path in entry.databases),
            }
            for entry in entries
        ],
        width="stretch",
        hide_index=True,
    )


def _download_table(
    label: str, rows: Sequence[Mapping[str, Any]], filename: str
) -> None:
    """Role: Compute download table for the app workflow."""
    st.download_button(
        label,
        data=records_to_csv(rows),
        file_name=filename,
        mime="text/csv",
        disabled=not rows,
    )


def _header(gene_count: int, database_count: int, master_path: Path) -> None:
    """Role: Compute header for the app workflow."""
    st.markdown(
        """
        <div class="hero">
          <div class="eyebrow">RARE DISEASE × ONCOLOGY</div>
          <h1>Gene Evidence Explorer</h1>
          <p>Interactive inspection of the immutable 2.1 molecular and therapy snapshots.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.caption(
        f"{gene_count} stored genes · {database_count} linked SQLite databases · "
        f"new executions → {_relative(master_path)}"
    )


def _gene_label(gene: str, entries: Mapping[str, Any]) -> str:
    """Role: Compute gene label for the app workflow."""
    if not gene:
        return "No stored genes"
    entry = entries[gene]
    therapy = " + therapy" if entry.has_therapy else ""
    return f"{gene} · {entry.analysis_runs} molecular run(s){therapy}"


def _relative(path: Path) -> str:
    """Role: Compute relative for the app workflow."""
    try:
        return str(path.resolve().relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def _count(payload: Mapping[str, Any] | None, key: str) -> int:
    """Role: Compute count for the app workflow."""
    if payload is None:
        return 0
    value = payload.get(key, 0)
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def _style() -> None:
    """Role: Compute style for the app workflow."""
    st.markdown(
        """
        <style>
        .stApp { background: linear-gradient(180deg, #f8fafc 0%, #fff7ed 100%); }
        .hero { padding: 1.4rem 1.6rem; border-radius: 18px; color: white;
                background: linear-gradient(115deg, #172033 0%, #7c2d12 58%, #ea580c 100%);
                box-shadow: 0 16px 40px rgba(23, 32, 51, .16); margin-bottom: .7rem; }
        .hero h1 { margin: .05rem 0 .25rem; font-size: 2.35rem; letter-spacing: -.03em; }
        .hero p { margin: 0; color: #ffedd5; }
        .eyebrow { font-size: .75rem; font-weight: 800; letter-spacing: .18em; color: #fdba74; }
        [data-testid="stMetric"] { background: rgba(255,255,255,.82); border: 1px solid #fed7aa;
                                    padding: .8rem; border-radius: 14px; }
        [data-testid="stSidebar"] { border-right: 1px solid #e2e8f0; }
        </style>
        """,
        unsafe_allow_html=True,
    )


if __name__ == "__main__":
    main()
