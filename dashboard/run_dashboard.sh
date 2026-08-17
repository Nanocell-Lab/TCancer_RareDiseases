#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="${DASHBOARD_VENV:-$HERE/.venv}"

if [[ ! -x "$VENV/bin/python" ]]; then
  python3 -m venv "$VENV"
fi

if ! "$VENV/bin/python" -c 'import streamlit, plotly, rare_oncology_dashboard' >/dev/null 2>&1; then
  "$VENV/bin/python" -m pip install --quiet -e "$HERE"
fi

exec "$VENV/bin/streamlit" run "$HERE/app.py" "$@"
