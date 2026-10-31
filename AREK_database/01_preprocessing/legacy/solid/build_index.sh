#!/usr/bin/env bash
set -euo pipefail
# Portable recovered source. Canonical resumable route: run_pipeline.py.
: "${RAW_DATA_ROOT:?set raw data root}" "${RESULT_ROOT:?set isolated output root}" "${BFAST_HOME:?set installed BFAST 0.5.6 prefix}" "${SOLID_ENV:?set dedicated colorspace environment prefix}" "${BFAST_DB_ROOT:?set BFAST hg38 CS database root}"
if [[ "${AREK_LEGACY_SOURCE_MODE:-}" != "isolated_research" ]]; then echo "Recovered historical entry; use canonical run_pipeline.py. To inspect native batching/diagnostics in a NEW output directory set AREK_LEGACY_SOURCE_MODE=isolated_research." >&2; exit 2; fi


# 避免旧环境动态库污染
unset LD_LIBRARY_PATH

BFAST_HOME="${BFAST_HOME}"
export PATH="${BFAST_HOME}/bin:/usr/bin:/bin:${PATH}"

HG38_FASTA="${HG38_FASTA:-${HG38_FASTA}}"

DB_ROOT="${BFAST_DB_ROOT}"
OUT_DIR="${DB_ROOT}/index"
TMP_DIR="${DB_ROOT}/tmp"
LOG_DIR="${DB_ROOT}/logs"

ID="hg38"

mkdir -p "${OUT_DIR}" "${TMP_DIR}" "${LOG_DIR}"

echo "[INFO] BFAST_HOME=${BFAST_HOME}"
echo "[INFO] HG38_FASTA=${HG38_FASTA}"
echo "[INFO] DB_ROOT=${DB_ROOT}"
echo "[INFO] OUT_DIR=${OUT_DIR}"
echo "[INFO] TMP_DIR=${TMP_DIR}"

if [[ ! -s "${HG38_FASTA}" ]]; then
  echo "[ERROR] HG38_FASTA not found or empty: ${HG38_FASTA}" >&2
  exit 1
fi

for cmd in bpreprocess bmatches balign bpostprocess; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "[ERROR] Missing command: ${cmd}" >&2
    exit 1
  fi
done

cat > "${OUT_DIR}/layouts.txt" <<'EOF'
10 111111111111111111111111
EOF

seq 0 100 | paste -sd ' ' - > "${OUT_DIR}/offsets.0-100.txt"

cat > "${OUT_DIR}/scoring.1.txt" <<'EOF'
-400
-70
100
-150
0
-125
EOF

echo "[INFO] Step 1: build NT reference brg"
if [[ ! -s "${OUT_DIR}/bfast.rg.file.${ID}.0.brg" ]]; then
  bpreprocess \
    -r "${HG38_FASTA}" \
    -a 0 \
    -A 0 \
    -o "${ID}" \
    -d "${OUT_DIR}/" \
    -T "${TMP_DIR}/" \
    > "${LOG_DIR}/01_bpreprocess_nt.log" 2>&1
else
  echo "[SKIP] NT brg exists"
fi

echo "[INFO] Step 2: build color-space reference brg"
if [[ ! -s "${OUT_DIR}/bfast.rg.file.${ID}.1.brg" ]]; then
  bpreprocess \
    -r "${HG38_FASTA}" \
    -a 0 \
    -A 1 \
    -o "${ID}" \
    -d "${OUT_DIR}/" \
    -T "${TMP_DIR}/" \
    > "${LOG_DIR}/02_bpreprocess_cs.log" 2>&1
else
  echo "[SKIP] CS brg exists"
fi

echo "[INFO] Step 3: build color-space BFAST index"
if ! ls "${OUT_DIR}"/*"${ID}"*.bif >/dev/null 2>&1; then
  bpreprocess \
    -r "${OUT_DIR}/bfast.rg.file.${ID}.1.brg" \
    -a 1 \
    -A 1 \
    -i "${OUT_DIR}/layouts.txt" \
    -o "${ID}" \
    -d "${OUT_DIR}/" \
    -T "${TMP_DIR}/" \
    > "${LOG_DIR}/03_bpreprocess_index.log" 2>&1
else
  echo "[SKIP] BIF index exists"
fi

ls -1 "${OUT_DIR}"/*"${ID}"*.bif > "${OUT_DIR}/main.indexes.${ID}.txt"
ls -1 "${OUT_DIR}"/*"${ID}"*.bif > "${OUT_DIR}/secondary.indexes.${ID}.txt"

echo "[INFO] Done. BFAST hg38 color-space index:"
ls -lh "${OUT_DIR}" | tee "${LOG_DIR}/final_index_files.txt"
