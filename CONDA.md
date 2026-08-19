# Conda environment

The repository provides one active environment for rare-oncology 2.1 and the
interactive dashboard. Python 3.11 is used because it satisfies every declared
`requires-python` constraint and is already exercised by the 2.0 and 2.1 test
environments.

## Create and activate

Run these commands from the repository root so the editable paths in
`environment.yml` resolve correctly:

```bash
conda env create --file environment.yml
conda activate rare-oncology-2.1
```

To update an existing environment after editing the file:

```bash
conda env update --file environment.yml --prune
```

## Validate

```bash
python -m unittest discover -s 2.1/tests -v
python -m unittest discover -s dashboard/tests -v
python -m rare_oncology.cli --help
streamlit run dashboard/app.py
```

The dashboard will use the active environment's `rare_oncology` module when a
dedicated `2.1/.venv` executable is absent. To force a particular CLI, set
`RARE_ONCOLOGY_CLI` to its absolute path.

## Why 2.0 is not installed

Both `2.0` and `2.1` declare the distribution name `rare-oncology` and expose
the import package `rare_oncology`. Installing both into one environment would
make import and console-script resolution order-dependent. Version 2.1 is the
active producer and supersedes 2.0. The historical 2.0 tests can still be run
without installing it:

```bash
PYTHONPATH=2.0/src python -m unittest discover -s 2.0/tests -v
```

## Dependency roles

- The core 2.1 data pipeline uses the Python standard library.
- Geneview, Matplotlib, pandas, SciPy and seaborn support protein-track plots.
- NumPy, Plotly and Streamlit support the interactive dashboard and structural
  comparisons.
- `setuptools`, `wheel` and `pip` install both repository packages in editable
  mode.

API credentials such as `ONCOKB_TOKEN`, `NCBI_API_KEY` and `NCBI_EMAIL` are not
stored in this file. Export them only in the shell when the corresponding source
is required.
