#!/usr/bin/env bash
set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GENE_LIST="${GENE_LIST:-${ROOT_DIR}/gene_lists/receptor_tyrosine_kinases.txt}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${ROOT_DIR}/results/rtk_maximum}"
DATABASE="${DATABASE:-${OUTPUT_ROOT}/rare_oncology.sqlite3}"
MAX_LOLLIPOP_RESIDUES="${MAX_LOLLIPOP_RESIDUES:-${MAX_PLOT_MUTATIONS:-80}}"
MAX_LOLLIPOP_LABELS="${MAX_LOLLIPOP_LABELS:-15}"
MAX_HEATMAP_MUTATIONS="${MAX_HEATMAP_MUTATIONS:-30}"
MAX_RARE_DISEASE_MUTATIONS="${MAX_RARE_DISEASE_MUTATIONS:-30}"
MAX_THERAPY_RECORDS="${MAX_THERAPY_RECORDS:-1000}"
STRING_PARTNERS="${STRING_PARTNERS:-20}"
STRING_SCORE="${STRING_SCORE:-700}"
STRING_NETWORK_TYPE="${STRING_NETWORK_TYPE:-functional}"
SLEEP_SECONDS="${SLEEP_SECONDS:-1}"

export MAX_LOLLIPOP_RESIDUES
export MAX_LOLLIPOP_LABELS
export MAX_HEATMAP_MUTATIONS
export MAX_RARE_DISEASE_MUTATIONS


if [[ -x "${ROOT_DIR}/.venv/bin/rare-oncology" ]]; then
  CLI="${ROOT_DIR}/.venv/bin/rare-oncology"
elif command -v rare-oncology >/dev/null 2>&1; then
  CLI="$(command -v rare-oncology)"
else
  echo "No se encontró rare-oncology. Instale 2.1 o cree 2.1/.venv." >&2
  exit 2
fi

if [[ ! -f "${GENE_LIST}" ]]; then
  echo "No existe la lista: ${GENE_LIST}" >&2
  exit 2
fi
if [[ -z "${NCBI_EMAIL:-}" ]]; then
  echo "Defina NCBI_EMAIL para las consultas E-utilities." >&2
  exit 2
fi
if [[ -z "${ONCOKB_TOKEN:-}" ]]; then
  echo "Defina ONCOKB_TOKEN para incluir OncoKB licenciado." >&2
  exit 2
fi

mkdir -p "${OUTPUT_ROOT}/logs"

common_args=(
  --program ALL
  --min-patients 1
  --classification ALL
  --source alphamissense
  --source gnomad
  --source panelapp
  --gnomad-max-variants 25
  --therapy-source civic
  --therapy-source oncokb
  --therapy-source clinicaltrials
  --therapy-source opentargets
  --max-therapy-records "${MAX_THERAPY_RECORDS}"
  --plots
  --plot-format png
  --plot-format svg
  --max-lollipop-residues "${MAX_LOLLIPOP_RESIDUES}"
  --max-lollipop-labels "${MAX_LOLLIPOP_LABELS}"
  --max-heatmap-mutations "${MAX_HEATMAP_MUTATIONS}"
  --max-rare-disease-mutations "${MAX_RARE_DISEASE_MUTATIONS}"
  --string-network
  --string-partners "${STRING_PARTNERS}"
  --string-score "${STRING_SCORE}"
  --string-network-type "${STRING_NETWORK_TYPE}"
  --database "${DATABASE}"
)

if [[ -n "${ORPHADATA_FILE:-}" ]]; then
  if [[ ! -f "${ORPHADATA_FILE}" ]]; then
    echo "ORPHADATA_FILE no existe: ${ORPHADATA_FILE}" >&2
    exit 2
  fi
  common_args+=(--orphadata-file "${ORPHADATA_FILE}")
fi
if [[ -n "${FDA_FILE:-}" ]]; then
  if [[ ! -f "${FDA_FILE}" ]]; then
    echo "FDA_FILE no existe: ${FDA_FILE}" >&2
    exit 2
  fi
  common_args+=(--regulatory-file "fda=${FDA_FILE}")
fi
if [[ -n "${EMA_FILE:-}" ]]; then
  if [[ ! -f "${EMA_FILE}" ]]; then
    echo "EMA_FILE no existe: ${EMA_FILE}" >&2
    exit 2
  fi
  common_args+=(--regulatory-file "ema=${EMA_FILE}")
fi
if [[ "${REFRESH:-0}" == "1" ]]; then
  common_args+=(--refresh)
fi

failures=()
while IFS= read -r gene || [[ -n "${gene}" ]]; do
  [[ -z "${gene}" ]] && continue
  if [[ ! "${gene}" =~ ^[A-Z][A-Z0-9-]*$ ]]; then
    echo "Símbolo inválido omitido: ${gene}" >&2
    failures+=("${gene}")
    continue
  fi

  gene_output="${OUTPUT_ROOT}/${gene}"
  mkdir -p "${gene_output}"
  echo "=== ${gene} ==="
  if "${CLI}" "${gene}" \
      "${common_args[@]}" \
      --output-dir "${gene_output}" \
      2>&1 | tee "${OUTPUT_ROOT}/logs/${gene}.log"; then
    echo "${gene}: completado"
  else
    echo "${gene}: falló; revise logs/${gene}.log" >&2
    failures+=("${gene}")
  fi
  sleep "${SLEEP_SECONDS}"
done < "${GENE_LIST}"

if (( ${#failures[@]} > 0 )); then
  printf 'Genes fallidos (%d): %s\n' "${#failures[@]}" "${failures[*]}" >&2
  exit 1
fi

echo "Completados todos los RTK. Resultados: ${OUTPUT_ROOT}"
