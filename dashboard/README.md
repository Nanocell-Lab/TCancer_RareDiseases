# Rare × Oncology Gene Dashboard

Independent Streamlit dashboard for the analyses stored by `2.1`.

A fresh clone opens with compact read-only examples for **EGFR, TGFBR1, ACVR1
and PDGFRB** from `../2.1/bootstrap/rare_oncology_demo.sqlite3`. New analyses
continue to use the writable master under `../2.1/results`; the seed is never
modified.

## What it does

- Opens on a read-only **Global overview** that selects the newest molecular
  snapshot per gene, links only therapy snapshots with the same stored molecular
  signature, and summarizes genes, protein families, actionability layers and
  rare-disease/cancer-project convergence without any web request.
- Filters the global atlas by curated receptor groups from `../2.1/gene_lists`,
  exports its gene/family projections, and drills directly into the selected
  gene's immutable evidence record.
- Discovers every compatible `*.sqlite3` under `../2.1/results` and reads it in
  SQLite read-only mode.
- Also discovers the tracked clone bootstrap and automatically adds its
  normalized FDA/EMA tables to newly requested therapy runs.
- Lists the genes already present and opens the newest molecular and therapy
  snapshots without calling the web.
- If an HGNC symbol is missing, **Load / run** calls the existing
  `../2.1/.venv/bin/rare-oncology` CLI. The 2.1 cache is checked first and any
  new execution is written to `../2.1/results/dashboard_master.sqlite3`.
- Offers interactive mutation overlap, protein lollipop, recurrence, oncology
  heatmap, top-30 rare-disease heatmap, AlphaMissense, three-layer therapy and
  human STRING partner-network views.
- Reads exact STRING snapshots from `string_network_cache`. A dedicated
  **Fetch/Refresh STRING** action runs only the STRING connector, so it does not
  repeat ClinVar, GDC or therapy requests.
- Generates validated mutant-sequence inputs for AlphaFold 3 and Boltz-2 from
  stored missense variants. Canonical sequences are reused from the 2.1 SQLite
  cache before optional UniProt retrieval.
- Compares uploaded wild-type and variant PDB files using a global C-alpha
  Kabsch alignment, interactive 3D overlay, RMSD and per-residue displacement.
- Preserves the raw 2.1 JSON, source statuses, request parameters and database
  location for provenance.

The existing RTK, TGF-β receptor and kinase-linked cytokine receptor databases
remain unchanged. The dashboard federates them rather than copying hundreds of
megabytes into another file.

## Start

```bash
cd dashboard
./run_dashboard.sh
```

The launcher creates `dashboard/.venv`, installs this package and starts
Streamlit. To pass Streamlit flags:

```bash
./run_dashboard.sh --server.port 8502
```

Manual installation is also supported:

```bash
cd dashboard
python3 -m venv .venv
.venv/bin/pip install -e .
.venv/bin/streamlit run app.py
```

## Optional environment variables

```bash
export NCBI_EMAIL="researcher@example.org"
export NCBI_API_KEY="..."       # optional NCBI rate-limit increase
export ONCOKB_TOKEN="..."       # licensed source; never shown in the UI/command
export RARE_ONCOLOGY_CLI="/absolute/path/to/rare-oncology"  # optional override
```

Do not commit tokens. The dashboard inherits them only in the child process;
they are not added to subprocess arguments or SQLite cache parameters.

## Cache-first behavior

1. All selected 2.1 SQLite files are queried for the gene.
2. If a stored snapshot exists, it is rendered immediately and no web call is
   made.
3. If it is absent, the user can run the full 2.1 analysis (`program=ALL`,
   `classification=ALL`) with AlphaMissense, gnomAD, PanelApp and a STRING
   partner network.
4. The 2.1 CLI performs its own exact-parameter cache lookup, calls only the
   missing web sources, and writes a new immutable snapshot to the master DB.
5. The dashboard reruns, discovers the master DB and renders the new snapshot.

For STRING specifically, the sidebar exposes the minimum score (0–1000), number
of added partners (0–100) and functional/physical network type. The partner tab
contains a bounded interactive network, direct-partner ranking, evidence-channel
heatmap and CSV/JSON downloads. Every direct gene–partner edge is retained in the
plot; only secondary partner–partner edges are capped to reduce visual clutter.

The same focused operation is available without Streamlit:

```bash
cd ../2.1
PYTHONPATH=src .venv/bin/python -m rare_oncology.string_cli TGFBR1 \
  --database results/dashboard_master.sqlite3 \
  --output-dir results/dashboard_runs/TGFBR1/string \
  --string-score 700 --string-partners 20 --string-network-type functional
```

Therapy results remain separated into variant-specific evidence, regulatory
therapy, and trial/target discovery. Nothing in this dashboard is a clinical
recommendation.

## Structural modelling

Open the **Structures** tab for either workflow:

1. Select one or more simple missense variants and optionally upload a canonical
   FASTA. The dashboard downloads one AlphaFold 3 JSON (dialect v4) and one
   Boltz-2 YAML per validated mutant, plus a manifest. The ZIP button bundles
   every generated input. Boltz-2 uses YAML as its native input format rather
   than AlphaFold's JSON schema.
2. Upload a wild-type PDB and a variant PDB. Only first-model `ATOM` C-alpha
   records are compared. Residues are paired by chain, residue number and
   insertion code, then the variant is superposed using Kabsch alignment.
   Results can be downloaded as JSON and residue-level CSV.

PDB uploads are parsed locally by the dashboard and are not sent to an external
service. RMSD is descriptive rather than proof of functional gain/loss; residue
numbering, missing regions, oligomer state, confidence and modelling protocol
must be comparable.

## Tests

```bash
cd dashboard
.venv/bin/python -m unittest discover -s tests -v
```

The tests use temporary SQLite fixtures and never access the network.

STRING data are association evidence, not proof of direct physical binding,
causality or therapeutic relevance. See the official [STRING API
documentation](https://version-12-0.string-db.org/help/api/).
