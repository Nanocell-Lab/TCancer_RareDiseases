# Clone bootstrap

This directory makes a fresh Git clone useful before any external API request.
It contains:

- `rare_oncology_demo.sqlite3`: compact, read-only molecular, therapy and STRING
  snapshots for **EGFR, TGFBR1, ACVR1 and PDGFRB**.
- `regulatory/fda_biomarkers.csv`: normalized FDA pharmacogenomic-label table.
- `regulatory/ema_labels.csv`: normalized EMA therapeutic-indication table.
- `manifest.json`: SHA-256 integrity, source execution identifiers and compaction
  policy.

The dashboard discovers the seed automatically in addition to writable
databases under `2.1/results`. New executions are always written to
`2.1/results/dashboard_master.sqlite3`; the tracked seed is never modified.

The molecular seed retains all shared mutations and all records needed for the
lollipop and heatmaps. Detailed ClinVar/GDC evidence rows are omitted to keep
the bundle below 10 MiB. AlphaMissense rows are deduplicated to protein changes
shown in the mutation plot, while gnomAD and PanelApp evidence remains. The
stored caveat identifies every compact snapshot.

Therapy layers remain semantically separate. The seed is for research and UI
demonstration, not clinical decision-making.

## Rebuild

After completing the full local analyses and regenerating FDA/EMA tables:

```bash
cd 2.1
python scripts/build_bootstrap_bundle.py
```

The builder selects the newest molecular analysis for each seed gene, requires
an exactly compatible therapy signature, copies the newest STRING network,
removes machine-local paths and enforces the 10 MiB limit.
