#!/usr/bin/env bash
set -Eeuo pipefail

###############################################################################
# Shotgun metagenome taxonomy profiling pipeline
# Default project: PRJNA613834
#
# Workflow:
#   1. Recursively scan and strictly validate paired-end FASTQ files
#   2. Sample read-length QC
#   3. fastp quality control and adapter trimming
#   4. Bowtie2 hg38 host depletion
#   5. Kraken2 taxonomic classification
#   6. Bracken abundance estimation at S/G/F/P levels
#   7. Merge sample outputs into feature/taxonomy tables
#   8. Optionally create BIOM and QIIME2 artifacts
#
# Main improvements:
#   - Stage-specific .done flags
#   - Atomic .part outputs
#   - FORCE=1 truly removes prior stage outputs
#   - Safe restart after interruption
#   - Faster fastp defaults for large shotgun FASTQ files
###############################################################################

############################
# User parameters
############################

PROJECT="${PROJECT:-PRJNA613834}"
RAW_ROOT="${RAW_ROOT:?Set RAW_ROOT to the original raw-read directory}"
OUT_ROOT="${OUT_ROOT:?Set OUT_ROOT to an isolated output directory}"
MISS_ROOT="${MISS_ROOT:?Set MISS_ROOT to an isolated output directory}"

KRAKEN_DB="${KRAKEN_DB:?Set KRAKEN_DB to the original Kraken2/Bracken database}"
HOST_BT2_DB="${HOST_BT2_DB:?Set HOST_BT2_DB to the original hg38 index prefix}"

SKIP_HOST_REMOVAL="${SKIP_HOST_REMOVAL:-0}"

# Keep the per-sample pipeline serial. Individual tools use these thread counts.
THREADS="${THREADS:-16}"
FASTP_THREADS="${FASTP_THREADS:-8}"
BOWTIE2_THREADS="${BOWTIE2_THREADS:-${THREADS}}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-8}"
PIGZ_THREADS="${PIGZ_THREADS:-8}"
KRAKEN_THREADS="${KRAKEN_THREADS:-${THREADS}}"

READ_LEN="${READ_LEN:-150}"
BRACKEN_THRESHOLD="${BRACKEN_THRESHOLD:-10}"
MIN_LENGTH="${MIN_LENGTH:-50}"

# fastp controls
FASTP_COMPRESSION="${FASTP_COMPRESSION:-1}"
FASTP_CORRECTION="${FASTP_CORRECTION:-0}"
FASTP_DISABLE_DUP_EVAL="${FASTP_DISABLE_DUP_EVAL:-1}"

# Optional raw gzip validation before fastp. This reads the entire input pair.
VERIFY_RAW_GZIP="${VERIFY_RAW_GZIP:-0}"

# 1: keep Kraken2 per-read output; 0: discard it.
KEEP_KRAKEN_OUTPUT="${KEEP_KRAKEN_OUTPUT:-0}"

# 1: remove clean/non-human FASTQ and temporary BAM after all Bracken stages.
CLEAN_INTERMEDIATE="${CLEAN_INTERMEDIATE:-1}"

# 1: archive an old QIIME2/16S-style project directory automatically.
ARCHIVE_OLD_16S="${ARCHIVE_OLD_16S:-0}"

# 1: delete all per-sample outputs and rerun all stages.
FORCE="${FORCE:-0}"

LENGTH_QC_READS="${LENGTH_QC_READS:-2000}"
TMP_ROOT="${TMP_ROOT:-${OUT_ROOT}/tmp}"
SAMTOOLS_SORT_MEM="${SAMTOOLS_SORT_MEM:-512M}"

BRACKEN_LEVELS=("S" "G" "F" "P")

# Existing standalone Kraken2 installation.


############################
# Helpers
############################

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

ACTIVE_SAMPLE=""
ACTIVE_STAGE=""
STATUS_DIR=""

mark_time() {
  local path="$1"
  mkdir -p "$(dirname "${path}")"
  date -Iseconds > "${path}"
}

fail() {
  local message="$*"

  if [[ -n "${ACTIVE_SAMPLE:-}" && -n "${STATUS_DIR:-}" ]]; then
    mark_time "${STATUS_DIR}/${ACTIVE_SAMPLE}.failed" 2>/dev/null || true
    rm -f "${STATUS_DIR}/${ACTIVE_SAMPLE}.running" 2>/dev/null || true

    if [[ -n "${ACTIVE_STAGE:-}" ]]; then
      mark_time "${STATUS_DIR}/${ACTIVE_SAMPLE}.${ACTIVE_STAGE}.failed" 2>/dev/null || true
      rm -f "${STATUS_DIR}/${ACTIVE_SAMPLE}.${ACTIVE_STAGE}.running" 2>/dev/null || true
    fi
  fi

  printf '[%s] [ERROR] %s\n' "$(date '+%F %T')" "${message}" >&2
  exit 1
}

on_error() {
  local exit_code=$?
  local line_no="${1:-unknown}"
  local command_text="${2:-unknown}"

  if [[ -n "${ACTIVE_SAMPLE:-}" && -n "${STATUS_DIR:-}" ]]; then
    mark_time "${STATUS_DIR}/${ACTIVE_SAMPLE}.failed" 2>/dev/null || true
    rm -f "${STATUS_DIR}/${ACTIVE_SAMPLE}.running" 2>/dev/null || true

    if [[ -n "${ACTIVE_STAGE:-}" ]]; then
      mark_time "${STATUS_DIR}/${ACTIVE_SAMPLE}.${ACTIVE_STAGE}.failed" 2>/dev/null || true
      rm -f "${STATUS_DIR}/${ACTIVE_SAMPLE}.${ACTIVE_STAGE}.running" 2>/dev/null || true
    fi
  fi

  printf '[%s] [ERROR] Pipeline failed at line %s, exit=%s, command=%s\n' \
    "$(date '+%F %T')" "${line_no}" "${exit_code}" "${command_text}" >&2
  exit "${exit_code}"
}

on_signal() {
  local signal_name="$1"

  if [[ -n "${ACTIVE_SAMPLE:-}" && -n "${STATUS_DIR:-}" ]]; then
    mark_time "${STATUS_DIR}/${ACTIVE_SAMPLE}.failed" 2>/dev/null || true
    rm -f "${STATUS_DIR}/${ACTIVE_SAMPLE}.running" 2>/dev/null || true

    if [[ -n "${ACTIVE_STAGE:-}" ]]; then
      mark_time "${STATUS_DIR}/${ACTIVE_SAMPLE}.${ACTIVE_STAGE}.failed" 2>/dev/null || true
      rm -f "${STATUS_DIR}/${ACTIVE_SAMPLE}.${ACTIVE_STAGE}.running" 2>/dev/null || true
    fi
  fi

  printf '[%s] [ERROR] Received signal %s\n' "$(date '+%F %T')" "${signal_name}" >&2
  exit 130
}

trap 'on_error "${LINENO}" "${BASH_COMMAND}"' ERR
trap 'on_signal INT' INT
trap 'on_signal TERM' TERM

require_bool() {
  local name="$1"
  local value="$2"
  case "${value}" in
    0|1) ;;
    *) fail "${name} must be 0 or 1; observed: ${value}" ;;
  esac
}

require_positive_integer() {
  local name="$1"
  local value="$2"
  [[ "${value}" =~ ^[1-9][0-9]*$ ]] || \
    fail "${name} must be a positive integer; observed: ${value}"
}

gzip_ok() {
  local path="$1"
  [[ -s "${path}" ]] && gzip -t "${path}" >/dev/null 2>&1
}

pair_ok() {
  gzip_ok "$1" && gzip_ok "$2"
}

stage_start() {
  local sample="$1"
  local stage="$2"
  ACTIVE_STAGE="${stage}"
  rm -f \
    "${STATUS_DIR}/${sample}.${stage}.done" \
    "${STATUS_DIR}/${sample}.${stage}.failed"
  mark_time "${STATUS_DIR}/${sample}.${stage}.running"
}

stage_finish() {
  local sample="$1"
  local stage="$2"
  mark_time "${STATUS_DIR}/${sample}.${stage}.done"
  rm -f \
    "${STATUS_DIR}/${sample}.${stage}.running" \
    "${STATUS_DIR}/${sample}.${stage}.failed"
  ACTIVE_STAGE=""
}

require_bool "SKIP_HOST_REMOVAL" "${SKIP_HOST_REMOVAL}"
require_bool "FASTP_CORRECTION" "${FASTP_CORRECTION}"
require_bool "FASTP_DISABLE_DUP_EVAL" "${FASTP_DISABLE_DUP_EVAL}"
require_bool "VERIFY_RAW_GZIP" "${VERIFY_RAW_GZIP}"
require_bool "KEEP_KRAKEN_OUTPUT" "${KEEP_KRAKEN_OUTPUT}"
require_bool "CLEAN_INTERMEDIATE" "${CLEAN_INTERMEDIATE}"
require_bool "ARCHIVE_OLD_16S" "${ARCHIVE_OLD_16S}"
require_bool "FORCE" "${FORCE}"

require_positive_integer "THREADS" "${THREADS}"
require_positive_integer "FASTP_THREADS" "${FASTP_THREADS}"
require_positive_integer "BOWTIE2_THREADS" "${BOWTIE2_THREADS}"
require_positive_integer "SAMTOOLS_THREADS" "${SAMTOOLS_THREADS}"
require_positive_integer "PIGZ_THREADS" "${PIGZ_THREADS}"
require_positive_integer "KRAKEN_THREADS" "${KRAKEN_THREADS}"
require_positive_integer "READ_LEN" "${READ_LEN}"
require_positive_integer "BRACKEN_THRESHOLD" "${BRACKEN_THRESHOLD}"
require_positive_integer "MIN_LENGTH" "${MIN_LENGTH}"
require_positive_integer "LENGTH_QC_READS" "${LENGTH_QC_READS}"

############################
# Protect against old 16S output
############################

if [[ -d "${OUT_ROOT}" ]]; then
  OLD_16S_MARKER=0

  if [[ -d "${OUT_ROOT}/01_batches" ]]; then
    OLD_16S_MARKER=1
  fi

  if find "${OUT_ROOT}" -maxdepth 3 -type f \
      \( -name '*demux.qza' -o -name 'single-end-demux.qza' -o -name 'paired-end-demux.qza' \) \
      -print -quit 2>/dev/null | grep -q .; then
    OLD_16S_MARKER=1
  fi

  if [[ "${OLD_16S_MARKER}" == "1" && ! -d "${OUT_ROOT}/03_kraken2_bracken" ]]; then
    if [[ "${ARCHIVE_OLD_16S}" == "1" ]]; then
      ARCHIVE_PATH="${OUT_ROOT}.old_16s_$(date +%Y%m%d_%H%M%S)"
      log "Old QIIME2/16S output detected; archiving to ${ARCHIVE_PATH}"
      mv "${OUT_ROOT}" "${ARCHIVE_PATH}"
    else
      fail "Old QIIME2/16S output detected in ${OUT_ROOT}. Set ARCHIVE_OLD_16S=1 or move it manually."
    fi
  fi
fi

############################
# Output directories
############################

MANIFEST_DIR="${OUT_ROOT}/00_manifest"
IMPORT_DIR="${OUT_ROOT}/01_import"
RAW_LINK_DIR="${IMPORT_DIR}/raw_fastq_links"

QC_HOST_DIR="${OUT_ROOT}/02_qc_host_removed"
CLEAN_DIR="${QC_HOST_DIR}/clean_fastq"
NONHUMAN_DIR="${QC_HOST_DIR}/nonhuman_fastq"
HOST_BAM_DIR="${QC_HOST_DIR}/bam"
QC_REPORT_DIR="${QC_HOST_DIR}/reports"

PROFILE_DIR="${OUT_ROOT}/03_kraken2_bracken"
KRAKEN_REPORT_DIR="${PROFILE_DIR}/kraken2_reports"
KRAKEN_OUTPUT_DIR="${PROFILE_DIR}/kraken2_outputs"
BRACKEN_DIR="${PROFILE_DIR}/bracken"

EXPORT_DIR="${OUT_ROOT}/04_export"
TAX_DIR="${OUT_ROOT}/05_taxonomy"
LOG_DIR="${OUT_ROOT}/logs"
STATUS_DIR="${OUT_ROOT}/status"
SCRIPT_DIR="${OUT_ROOT}/scripts"

mkdir -p \
  "${MANIFEST_DIR}" \
  "${RAW_LINK_DIR}" \
  "${CLEAN_DIR}" \
  "${NONHUMAN_DIR}" \
  "${HOST_BAM_DIR}" \
  "${QC_REPORT_DIR}" \
  "${KRAKEN_REPORT_DIR}" \
  "${KRAKEN_OUTPUT_DIR}" \
  "${BRACKEN_DIR}" \
  "${EXPORT_DIR}" \
  "${TAX_DIR}" \
  "${LOG_DIR}" \
  "${STATUS_DIR}" \
  "${SCRIPT_DIR}" \
  "${TMP_ROOT}" \
  "${MISS_ROOT}"

for LVL in "${BRACKEN_LEVELS[@]}"; do
  mkdir -p "${BRACKEN_DIR}/${LVL}"
done

MANIFEST="${MANIFEST_DIR}/manifest.tsv"
MANIFEST_PAIRED="${MANIFEST_DIR}/manifest_paired.tsv"
FASTQ_USED="${MANIFEST_DIR}/fastq_used.tsv"
ALL_CLASS="${MANIFEST_DIR}/all_fastq_classification.tsv"
READ_LAYOUT="${MANIFEST_DIR}/read_layout.tsv"
MISSING_FASTQ="${MANIFEST_DIR}/missing_or_ambiguous_fastq_pairs.tsv"
LENGTH_QC="${MANIFEST_DIR}/pair_mate_length_qc.tsv"
PROJECT_INFO="${MANIFEST_DIR}/project_run_parameters.tsv"
MISS_COPY="${MISS_ROOT}/${PROJECT}_metagenome_missing_or_ambiguous_fastq_pairs.tsv"

############################
# Dependencies and databases
############################

log "=================================================================="
log "Project: ${PROJECT}"
log "RAW_ROOT: ${RAW_ROOT}"
log "OUT_ROOT: ${OUT_ROOT}"
log "KRAKEN_DB: ${KRAKEN_DB}"
log "HOST_BT2_DB: ${HOST_BT2_DB}"
log "SKIP_HOST_REMOVAL: ${SKIP_HOST_REMOVAL}"
log "THREADS: ${THREADS}"
log "FASTP_THREADS: ${FASTP_THREADS}"
log "BOWTIE2_THREADS: ${BOWTIE2_THREADS}"
log "SAMTOOLS_THREADS: ${SAMTOOLS_THREADS}"
log "PIGZ_THREADS: ${PIGZ_THREADS}"
log "KRAKEN_THREADS: ${KRAKEN_THREADS}"
log "READ_LEN: ${READ_LEN}"
log "MIN_LENGTH: ${MIN_LENGTH}"
log "FASTP_COMPRESSION: ${FASTP_COMPRESSION}"
log "FASTP_CORRECTION: ${FASTP_CORRECTION}"
log "FASTP_DISABLE_DUP_EVAL: ${FASTP_DISABLE_DUP_EVAL}"
log "VERIFY_RAW_GZIP: ${VERIFY_RAW_GZIP}"
log "KEEP_KRAKEN_OUTPUT: ${KEEP_KRAKEN_OUTPUT}"
log "CLEAN_INTERMEDIATE: ${CLEAN_INTERMEDIATE}"
log "FORCE: ${FORCE}"
log "=================================================================="

[[ -d "${RAW_ROOT}" ]] || fail "RAW_ROOT does not exist: ${RAW_ROOT}"
[[ -d "${KRAKEN_DB}" ]] || fail "KRAKEN_DB does not exist: ${KRAKEN_DB}"

for f in hash.k2d opts.k2d taxo.k2d; do
  [[ -s "${KRAKEN_DB}/${f}" ]] || fail "Missing Kraken2 database file: ${KRAKEN_DB}/${f}"
done

[[ -s "${KRAKEN_DB}/database${READ_LEN}mers.kmer_distrib" ]] || \
  fail "Missing Bracken distribution: ${KRAKEN_DB}/database${READ_LEN}mers.kmer_distrib"

if [[ "${SKIP_HOST_REMOVAL}" == "0" ]]; then
  compgen -G "${HOST_BT2_DB}.*.bt2*" >/dev/null || \
    fail "Cannot find Bowtie2 index with prefix: ${HOST_BT2_DB}"
fi

REQUIRED_COMMANDS=(fastp samtools pigz kraken2 bracken python3 gzip)
if [[ "${SKIP_HOST_REMOVAL}" == "0" ]]; then
  REQUIRED_COMMANDS+=(bowtie2)
fi

for cmd in "${REQUIRED_COMMANDS[@]}"; do
  command -v "${cmd}" >/dev/null 2>&1 || fail "Required command not found: ${cmd}"
done

python3 - <<'PY'
import importlib.util
import sys

for module in ("pandas",):
    if importlib.util.find_spec(module) is None:
        print(f"[ERROR] Required Python module not found: {module}", file=sys.stderr)
        sys.exit(1)
PY

{
  printf 'parameter\tvalue\n'
  printf 'project\t%s\n' "${PROJECT}"
  printf 'raw_root\t%s\n' "${RAW_ROOT}"
  printf 'out_root\t%s\n' "${OUT_ROOT}"
  printf 'kraken_db\t%s\n' "${KRAKEN_DB}"
  printf 'host_bt2_db\t%s\n' "${HOST_BT2_DB}"
  printf 'skip_host_removal\t%s\n' "${SKIP_HOST_REMOVAL}"
  printf 'threads\t%s\n' "${THREADS}"
  printf 'fastp_threads\t%s\n' "${FASTP_THREADS}"
  printf 'bowtie2_threads\t%s\n' "${BOWTIE2_THREADS}"
  printf 'samtools_threads\t%s\n' "${SAMTOOLS_THREADS}"
  printf 'pigz_threads\t%s\n' "${PIGZ_THREADS}"
  printf 'kraken_threads\t%s\n' "${KRAKEN_THREADS}"
  printf 'read_len\t%s\n' "${READ_LEN}"
  printf 'min_length\t%s\n' "${MIN_LENGTH}"
  printf 'bracken_threshold\t%s\n' "${BRACKEN_THRESHOLD}"
  printf 'fastp_compression\t%s\n' "${FASTP_COMPRESSION}"
  printf 'fastp_correction\t%s\n' "${FASTP_CORRECTION}"
  printf 'fastp_disable_dup_eval\t%s\n' "${FASTP_DISABLE_DUP_EVAL}"
  printf 'verify_raw_gzip\t%s\n' "${VERIFY_RAW_GZIP}"
  printf 'keep_kraken_output\t%s\n' "${KEEP_KRAKEN_OUTPUT}"
  printf 'clean_intermediate\t%s\n' "${CLEAN_INTERMEDIATE}"
  printf 'force\t%s\n' "${FORCE}"
  printf 'run_start\t%s\n' "$(date -Iseconds)"
} > "${PROJECT_INFO}"

###############################################################################
# 00_manifest: scan and strictly validate PE FASTQ layout
###############################################################################

log "Scanning FASTQ files and validating paired-end layout..."

if python3 - <<PY
from pathlib import Path
import os
import re
import sys
import pandas as pd

raw_root = Path(r"${RAW_ROOT}")
manifest = Path(r"${MANIFEST}")
manifest_paired = Path(r"${MANIFEST_PAIRED}")
fastq_used = Path(r"${FASTQ_USED}")
all_class = Path(r"${ALL_CLASS}")
read_layout = Path(r"${READ_LAYOUT}")
missing = Path(r"${MISSING_FASTQ}")
raw_link_dir = Path(r"${RAW_LINK_DIR}")
raw_link_dir.mkdir(parents=True, exist_ok=True)

fastq_exts = [".fastq.gz", ".fq.gz", ".fastq", ".fq"]
patterns = [
    (re.compile(r"(.+)_R1_001$", re.I), "R1"),
    (re.compile(r"(.+)_R2_001$", re.I), "R2"),
    (re.compile(r"(.+)_R1$", re.I), "R1"),
    (re.compile(r"(.+)_R2$", re.I), "R2"),
    (re.compile(r"(.+)\\.R1$", re.I), "R1"),
    (re.compile(r"(.+)\\.R2$", re.I), "R2"),
    (re.compile(r"(.+)_1$", re.I), "R1"),
    (re.compile(r"(.+)_2$", re.I), "R2"),
    (re.compile(r"(.+)\\.1$", re.I), "R1"),
    (re.compile(r"(.+)\\.2$", re.I), "R2"),
    (re.compile(r"(.+)_f1$", re.I), "R1"),
    (re.compile(r"(.+)_r1$", re.I), "R1"),
    (re.compile(r"(.+)_f2$", re.I), "R2"),
    (re.compile(r"(.+)_r2$", re.I), "R2"),
]

def matched_ext(name: str) -> str:
    lower = name.lower()
    for ext in fastq_exts:
        if lower.endswith(ext):
            return name[-len(ext):]
    return ""

def strip_fastq_ext(name: str) -> str:
    ext = matched_ext(name)
    return name[:-len(ext)] if ext else name

def classify_fastq(path: Path):
    stem = strip_fastq_ext(path.name)
    for pat, mate in patterns:
        match = pat.fullmatch(stem)
        if match:
            return match.group(1), mate
    return stem, "UNKNOWN"

records = {}
all_rows = []

for path in sorted(raw_root.rglob("*")):
    if not path.is_file():
        continue
    if path.name.endswith(".part"):
        continue
    lower = path.name.lower()
    if not any(lower.endswith(ext) for ext in fastq_exts):
        continue
    if path.stat().st_size == 0:
        continue

    sample, mate = classify_fastq(path)
    records.setdefault(sample, {"R1": [], "R2": [], "UNKNOWN": []})
    records[sample][mate].append(str(path.resolve()))
    all_rows.append(
        {
            "sample_id_guess": sample,
            "mate": mate,
            "fastq_path": str(path.resolve()),
            "file_size": path.stat().st_size,
        }
    )

all_df = pd.DataFrame(
    all_rows,
    columns=["sample_id_guess", "mate", "fastq_path", "file_size"],
)
all_df.to_csv(all_class, sep="\t", index=False)

ok_rows = []
bad_rows = []
used_rows = []
layout_rows = []

for sample, data in sorted(records.items()):
    r1s = data["R1"]
    r2s = data["R2"]
    unknowns = data["UNKNOWN"]

    if len(r1s) == 1 and len(r2s) == 1 and len(unknowns) == 0:
        ok_rows.append(
            {
                "sample_id": sample,
                "layout": "PE",
                "r1": r1s[0],
                "r2": r2s[0],
                "status": "OK",
            }
        )
        used_rows.extend(
            [
                {"sample_id": sample, "mate": "R1", "fastq_path": r1s[0]},
                {"sample_id": sample, "mate": "R2", "fastq_path": r2s[0]},
            ]
        )
        layout_rows.append({"sample_id": sample, "layout": "PE"})

        for mate, source in (("R1", r1s[0]), ("R2", r2s[0])):
            source_path = Path(source)
            suffix = matched_ext(source_path.name)
            link = raw_link_dir / f"{sample}_{mate}{suffix}"

            if link.exists() or link.is_symlink():
                try:
                    if link.resolve() == source_path.resolve():
                        continue
                except FileNotFoundError:
                    pass
                link.unlink()

            os.symlink(source, link)
    else:
        bad_rows.append(
            {
                "sample_id": sample,
                "r1": ";".join(r1s),
                "r2": ";".join(r2s),
                "unknown": ";".join(unknowns),
                "status": (
                    f"R1_count={len(r1s)}|"
                    f"R2_count={len(r2s)}|"
                    f"UNKNOWN_count={len(unknowns)}"
                ),
            }
        )

ok_df = pd.DataFrame(
    ok_rows,
    columns=["sample_id", "layout", "r1", "r2", "status"],
)
bad_df = pd.DataFrame(
    bad_rows,
    columns=["sample_id", "r1", "r2", "unknown", "status"],
)

ok_df.to_csv(manifest, sep="\t", index=False)
ok_df.to_csv(manifest_paired, sep="\t", index=False)
pd.DataFrame(
    used_rows,
    columns=["sample_id", "mate", "fastq_path"],
).to_csv(fastq_used, sep="\t", index=False)
pd.DataFrame(
    layout_rows,
    columns=["sample_id", "layout"],
).to_csv(read_layout, sep="\t", index=False)
bad_df.to_csv(missing, sep="\t", index=False)

print(f"[INFO] FASTQ files detected: {len(all_df)}")
print(f"[INFO] Strict paired-end samples: {len(ok_df)}")
print(f"[INFO] Missing/ambiguous/unpaired entries: {len(bad_df)}")

if len(all_df) == 0:
    sys.exit(2)
if len(ok_df) == 0:
    sys.exit(3)
if len(bad_df) > 0:
    sys.exit(4)
PY
then
  SCAN_RC=0
else
  SCAN_RC=$?
fi

cp -f "${MISSING_FASTQ}" "${MISS_COPY}"

case "${SCAN_RC}" in
  0) ;;
  2) fail "No FASTQ files detected under ${RAW_ROOT}." ;;
  3) fail "No strict paired-end samples detected. See ${MISSING_FASTQ}" ;;
  4) fail "Missing, ambiguous, duplicated, or unpaired FASTQ entries detected. See ${MISSING_FASTQ}" ;;
  *) fail "FASTQ manifest construction failed with exit code ${SCAN_RC}." ;;
esac

###############################################################################
# Read-length QC
###############################################################################

log "Running read-length QC using up to ${LENGTH_QC_READS} reads per mate..."

python3 - <<PY
from pathlib import Path
import gzip
import statistics
import pandas as pd
import sys

manifest = pd.read_csv(r"${MANIFEST}", sep="\t")
out = Path(r"${LENGTH_QC}")
max_reads = int("${LENGTH_QC_READS}")
min_required = int("${MIN_LENGTH}")

def opener(path):
    path = str(path)
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")

def lengths(path):
    values = []
    with opener(path) as handle:
        for index, line in enumerate(handle):
            if index % 4 == 1:
                values.append(len(line.rstrip("\r\n")))
                if len(values) >= max_reads:
                    break
    return values

rows = []
bad = []

for row in manifest.itertuples(index=False):
    r1_len = lengths(row.r1)
    r2_len = lengths(row.r2)

    if not r1_len or not r2_len:
        bad.append((row.sample_id, "empty_or_invalid_fastq"))
        continue

    record = {
        "sample_id": str(row.sample_id),
        "r1_reads_checked": len(r1_len),
        "r1_min": min(r1_len),
        "r1_median": statistics.median(r1_len),
        "r1_max": max(r1_len),
        "r2_reads_checked": len(r2_len),
        "r2_min": min(r2_len),
        "r2_median": statistics.median(r2_len),
        "r2_max": max(r2_len),
    }
    rows.append(record)

    if record["r1_median"] < min_required or record["r2_median"] < min_required:
        bad.append((row.sample_id, f"median_below_{min_required}"))

pd.DataFrame(rows).to_csv(out, sep="\t", index=False)

if bad:
    print("[ERROR] Read-length QC failed:", file=sys.stderr)
    for sample, reason in bad[:20]:
        print(f"  {sample}\t{reason}", file=sys.stderr)
    print(f"[ERROR] Full table: {out}", file=sys.stderr)
    sys.exit(5)

medians = [x["r1_median"] for x in rows] + [x["r2_median"] for x in rows]
print(f"[INFO] Samples checked: {len(rows)}")
print(f"[INFO] Global median of mate medians: {statistics.median(medians):.1f}")
print(f"[INFO] Read-length QC table: {out}")
PY

###############################################################################
# Per-sample helpers
###############################################################################

sample_outputs_complete() {
  local sample="$1"
  local level

  [[ -s "${KRAKEN_REPORT_DIR}/${sample}.kreport" ]] || return 1

  for level in "${BRACKEN_LEVELS[@]}"; do
    [[ -s "${BRACKEN_DIR}/${level}/${sample}.bracken.tsv" ]] || return 1
  done

  return 0
}

reset_sample_outputs() {
  local sample="$1"
  local level

  log "FORCE=1: removing previous outputs for ${sample}"

  rm -f \
    "${STATUS_DIR}/${sample}.done" \
    "${STATUS_DIR}/${sample}.running" \
    "${STATUS_DIR}/${sample}.failed" \
    "${STATUS_DIR}/${sample}.fastp.done" \
    "${STATUS_DIR}/${sample}.fastp.running" \
    "${STATUS_DIR}/${sample}.fastp.failed" \
    "${STATUS_DIR}/${sample}.host.done" \
    "${STATUS_DIR}/${sample}.host.running" \
    "${STATUS_DIR}/${sample}.host.failed" \
    "${STATUS_DIR}/${sample}.kraken.done" \
    "${STATUS_DIR}/${sample}.kraken.running" \
    "${STATUS_DIR}/${sample}.kraken.failed" \
    "${STATUS_DIR}/${sample}.bracken.done" \
    "${STATUS_DIR}/${sample}.bracken.running" \
    "${STATUS_DIR}/${sample}.bracken.failed"

  rm -f \
    "${CLEAN_DIR}/${sample}_R1.clean.fastq.gz" \
    "${CLEAN_DIR}/${sample}_R2.clean.fastq.gz" \
    "${CLEAN_DIR}/.${sample}_R1.clean.part.fastq.gz" \
    "${CLEAN_DIR}/.${sample}_R2.clean.part.fastq.gz" \
    "${NONHUMAN_DIR}/${sample}_R1.nonhuman.fastq.gz" \
    "${NONHUMAN_DIR}/${sample}_R2.nonhuman.fastq.gz" \
    "${NONHUMAN_DIR}/.${sample}_R1.nonhuman.part.fastq.gz" \
    "${NONHUMAN_DIR}/.${sample}_R2.nonhuman.part.fastq.gz" \
    "${HOST_BAM_DIR}/${sample}.both_unmapped.name_sorted.bam" \
    "${HOST_BAM_DIR}/.${sample}.both_unmapped.name_sorted.part.bam" \
    "${QC_REPORT_DIR}/${sample}.fastp.html" \
    "${QC_REPORT_DIR}/${sample}.fastp.json" \
    "${QC_REPORT_DIR}/.${sample}.fastp.part.html" \
    "${QC_REPORT_DIR}/.${sample}.fastp.part.json" \
    "${KRAKEN_REPORT_DIR}/${sample}.kreport" \
    "${KRAKEN_REPORT_DIR}/.${sample}.part.kreport" \
    "${KRAKEN_OUTPUT_DIR}/${sample}.kraken.out" \
    "${KRAKEN_OUTPUT_DIR}/.${sample}.part.kraken.out"

  for level in "${BRACKEN_LEVELS[@]}"; do
    rm -f \
      "${BRACKEN_DIR}/${level}/${sample}.bracken.tsv" \
      "${BRACKEN_DIR}/${level}/.${sample}.part.bracken.tsv" \
      "${STATUS_DIR}/${sample}.bracken.${level}.done"
  done

  rm -rf "${TMP_ROOT:?}/${sample}.work"
}

###############################################################################
# Per-sample processing
###############################################################################

process_sample() {
  local SAMPLE_ID="$1"
  local R1="$2"
  local R2="$3"

  local SAMPLE_DONE="${STATUS_DIR}/${SAMPLE_ID}.done"
  local SAMPLE_RUNNING="${STATUS_DIR}/${SAMPLE_ID}.running"
  local SAMPLE_FAILED="${STATUS_DIR}/${SAMPLE_ID}.failed"

  local FASTP_DONE="${STATUS_DIR}/${SAMPLE_ID}.fastp.done"
  local HOST_DONE="${STATUS_DIR}/${SAMPLE_ID}.host.done"
  local KRAKEN_DONE="${STATUS_DIR}/${SAMPLE_ID}.kraken.done"
  local BRACKEN_DONE="${STATUS_DIR}/${SAMPLE_ID}.bracken.done"

  local CLEAN_R1="${CLEAN_DIR}/${SAMPLE_ID}_R1.clean.fastq.gz"
  local CLEAN_R2="${CLEAN_DIR}/${SAMPLE_ID}_R2.clean.fastq.gz"
  local CLEAN_R1_PART="${CLEAN_DIR}/.${SAMPLE_ID}_R1.clean.part.fastq.gz"
  local CLEAN_R2_PART="${CLEAN_DIR}/.${SAMPLE_ID}_R2.clean.part.fastq.gz"

  local FASTP_HTML="${QC_REPORT_DIR}/${SAMPLE_ID}.fastp.html"
  local FASTP_JSON="${QC_REPORT_DIR}/${SAMPLE_ID}.fastp.json"
  local FASTP_HTML_PART="${QC_REPORT_DIR}/.${SAMPLE_ID}.fastp.part.html"
  local FASTP_JSON_PART="${QC_REPORT_DIR}/.${SAMPLE_ID}.fastp.part.json"

  local NONHUMAN_R1="${NONHUMAN_DIR}/${SAMPLE_ID}_R1.nonhuman.fastq.gz"
  local NONHUMAN_R2="${NONHUMAN_DIR}/${SAMPLE_ID}_R2.nonhuman.fastq.gz"
  local NONHUMAN_R1_PART="${NONHUMAN_DIR}/.${SAMPLE_ID}_R1.nonhuman.part.fastq.gz"
  local NONHUMAN_R2_PART="${NONHUMAN_DIR}/.${SAMPLE_ID}_R2.nonhuman.part.fastq.gz"

  local HOST_BAM="${HOST_BAM_DIR}/${SAMPLE_ID}.both_unmapped.name_sorted.bam"
  local HOST_BAM_PART="${HOST_BAM_DIR}/.${SAMPLE_ID}.both_unmapped.name_sorted.part.bam"

  local KREPORT="${KRAKEN_REPORT_DIR}/${SAMPLE_ID}.kreport"
  local KREPORT_PART="${KRAKEN_REPORT_DIR}/.${SAMPLE_ID}.part.kreport"

  local KOUTPUT=""
  local KOUTPUT_PART=""
  local WORK_DIR="${TMP_ROOT}/${SAMPLE_ID}.work"

  if [[ "${KEEP_KRAKEN_OUTPUT}" == "1" ]]; then
    KOUTPUT="${KRAKEN_OUTPUT_DIR}/${SAMPLE_ID}.kraken.out"
    KOUTPUT_PART="${KRAKEN_OUTPUT_DIR}/.${SAMPLE_ID}.part.kraken.out"
  else
    KOUTPUT="/dev/null"
    KOUTPUT_PART="/dev/null"
  fi

  if [[ "${FORCE}" == "1" ]]; then
    reset_sample_outputs "${SAMPLE_ID}"
  elif [[ -s "${SAMPLE_DONE}" ]] && sample_outputs_complete "${SAMPLE_ID}"; then
    log "Skip completed sample: ${SAMPLE_ID}"
    return 0
  elif [[ -s "${SAMPLE_DONE}" ]]; then
    log "Sample .done flag exists but final outputs are incomplete; resuming ${SAMPLE_ID}"
    rm -f "${SAMPLE_DONE}"
  fi

  rm -f "${SAMPLE_FAILED}"
  mark_time "${SAMPLE_RUNNING}"
  mkdir -p "${WORK_DIR}"

  log "=================================================================="
  log "Processing sample: ${SAMPLE_ID}"
  log "R1: ${R1}"
  log "R2: ${R2}"
  log "=================================================================="

  ############################
  # Optional raw gzip validation
  ############################

  if [[ "${VERIFY_RAW_GZIP}" == "1" ]]; then
    log "Validating raw gzip files: ${SAMPLE_ID}"
    [[ "${R1,,}" == *.gz ]] && gzip -t "${R1}"
    [[ "${R2,,}" == *.gz ]] && gzip -t "${R2}"
  fi

  ############################
  # fastp PE
  ############################

  if [[ -s "${FASTP_DONE}" ]] && ! pair_ok "${CLEAN_R1}" "${CLEAN_R2}"; then
    log "fastp .done flag is inconsistent with outputs; rebuilding ${SAMPLE_ID}"
    rm -f "${FASTP_DONE}" "${CLEAN_R1}" "${CLEAN_R2}"
  fi

  if [[ ! -s "${FASTP_DONE}" ]]; then
    stage_start "${SAMPLE_ID}" "fastp"
    log "fastp PE: ${SAMPLE_ID}"

    rm -f \
      "${CLEAN_R1_PART}" "${CLEAN_R2_PART}" \
      "${FASTP_HTML_PART}" "${FASTP_JSON_PART}"

    FASTP_ARGS=(
      --thread "${FASTP_THREADS}"
      --in1 "${R1}"
      --in2 "${R2}"
      --out1 "${CLEAN_R1_PART}"
      --out2 "${CLEAN_R2_PART}"
      --detect_adapter_for_pe
      --length_required "${MIN_LENGTH}"
      --compression "${FASTP_COMPRESSION}"
      --html "${FASTP_HTML_PART}"
      --json "${FASTP_JSON_PART}"
    )

    if [[ "${FASTP_CORRECTION}" == "1" ]]; then
      FASTP_ARGS+=(--correction)
    fi

    if [[ "${FASTP_DISABLE_DUP_EVAL}" == "1" ]] && \
       fastp --help 2>&1 | grep -q -- '--dont_eval_duplication'; then
      FASTP_ARGS+=(--dont_eval_duplication)
    fi

    fastp "${FASTP_ARGS[@]}" \
      > "${LOG_DIR}/${SAMPLE_ID}.fastp.log" 2>&1

    pair_ok "${CLEAN_R1_PART}" "${CLEAN_R2_PART}" || \
      fail "fastp produced invalid gzip output for ${SAMPLE_ID}"

    [[ -s "${FASTP_HTML_PART}" ]] || fail "Missing fastp HTML for ${SAMPLE_ID}"
    [[ -s "${FASTP_JSON_PART}" ]] || fail "Missing fastp JSON for ${SAMPLE_ID}"

    mv -f "${CLEAN_R1_PART}" "${CLEAN_R1}"
    mv -f "${CLEAN_R2_PART}" "${CLEAN_R2}"
    mv -f "${FASTP_HTML_PART}" "${FASTP_HTML}"
    mv -f "${FASTP_JSON_PART}" "${FASTP_JSON}"

    stage_finish "${SAMPLE_ID}" "fastp"
  else
    log "Completed fastp output detected: ${SAMPLE_ID}"
  fi

  pair_ok "${CLEAN_R1}" "${CLEAN_R2}" || \
    fail "Invalid cleaned FASTQ pair: ${SAMPLE_ID}"

  if [[ "${AC_STOP_AFTER_QC:-0}" == "1" ]]; then
    log "fastp preprocessing completed; host depletion and profiling were not run."
    return 0
  fi

  ############################
  # Host removal PE
  ############################

  local KRAKEN_R1
  local KRAKEN_R2

  if [[ "${SKIP_HOST_REMOVAL}" == "0" ]]; then
    if [[ -s "${HOST_DONE}" ]] && ! pair_ok "${NONHUMAN_R1}" "${NONHUMAN_R2}"; then
      log "Host-removal .done flag is inconsistent with outputs; rebuilding ${SAMPLE_ID}"
      rm -f "${HOST_DONE}" "${NONHUMAN_R1}" "${NONHUMAN_R2}"
    fi

    if [[ ! -s "${HOST_DONE}" ]]; then
      stage_start "${SAMPLE_ID}" "host"
      log "Bowtie2 host removal PE: ${SAMPLE_ID}"

      rm -f \
        "${HOST_BAM_PART}" \
        "${NONHUMAN_R1_PART}" \
        "${NONHUMAN_R2_PART}"

      bowtie2 \
        -x "${HOST_BT2_DB}" \
        -1 "${CLEAN_R1}" \
        -2 "${CLEAN_R2}" \
        --very-sensitive \
        --threads "${BOWTIE2_THREADS}" \
        2> "${LOG_DIR}/${SAMPLE_ID}.bowtie2_host.log" \
      | samtools view \
          -@ "${SAMTOOLS_THREADS}" \
          -b \
          -f 12 \
          -F 256 \
          - \
      | samtools sort \
          -@ "${SAMTOOLS_THREADS}" \
          -m "${SAMTOOLS_SORT_MEM}" \
          -n \
          -T "${WORK_DIR}/${SAMPLE_ID}.sorttmp" \
          -o "${HOST_BAM_PART}" \
          -

      samtools quickcheck -v "${HOST_BAM_PART}"

      local FIFO_R1="${WORK_DIR}/${SAMPLE_ID}.R1.fifo"
      local FIFO_R2="${WORK_DIR}/${SAMPLE_ID}.R2.fifo"
      rm -f "${FIFO_R1}" "${FIFO_R2}"
      mkfifo "${FIFO_R1}" "${FIFO_R2}"

      pigz -p "${PIGZ_THREADS}" -1 -c < "${FIFO_R1}" > "${NONHUMAN_R1_PART}" &
      local PIGZ_PID1=$!

      pigz -p "${PIGZ_THREADS}" -1 -c < "${FIFO_R2}" > "${NONHUMAN_R2_PART}" &
      local PIGZ_PID2=$!

      local SAMTOOLS_FASTQ_RC=0
      local PIGZ_RC1=0
      local PIGZ_RC2=0

      if samtools fastq \
        -@ "${SAMTOOLS_THREADS}" \
        -1 "${FIFO_R1}" \
        -2 "${FIFO_R2}" \
        -0 /dev/null \
        -s /dev/null \
        -n \
        "${HOST_BAM_PART}" \
        > "${LOG_DIR}/${SAMPLE_ID}.samtools_fastq.log" 2>&1
      then
        SAMTOOLS_FASTQ_RC=0
      else
        SAMTOOLS_FASTQ_RC=$?
      fi

      if wait "${PIGZ_PID1}"; then
        PIGZ_RC1=0
      else
        PIGZ_RC1=$?
      fi

      if wait "${PIGZ_PID2}"; then
        PIGZ_RC2=0
      else
        PIGZ_RC2=$?
      fi

      rm -f "${FIFO_R1}" "${FIFO_R2}"

      [[ "${SAMTOOLS_FASTQ_RC}" == "0" ]] || \
        fail "samtools fastq failed for ${SAMPLE_ID}, exit=${SAMTOOLS_FASTQ_RC}"
      [[ "${PIGZ_RC1}" == "0" && "${PIGZ_RC2}" == "0" ]] || \
        fail "pigz failed for ${SAMPLE_ID}: R1=${PIGZ_RC1}, R2=${PIGZ_RC2}"

      pair_ok "${NONHUMAN_R1_PART}" "${NONHUMAN_R2_PART}" || \
        fail "Host removal produced invalid non-human FASTQ for ${SAMPLE_ID}"

      mv -f "${HOST_BAM_PART}" "${HOST_BAM}"
      mv -f "${NONHUMAN_R1_PART}" "${NONHUMAN_R1}"
      mv -f "${NONHUMAN_R2_PART}" "${NONHUMAN_R2}"

      stage_finish "${SAMPLE_ID}" "host"
    else
      log "Completed host-removal output detected: ${SAMPLE_ID}"
    fi

    pair_ok "${NONHUMAN_R1}" "${NONHUMAN_R2}" || \
      fail "Invalid non-human FASTQ pair: ${SAMPLE_ID}"

    KRAKEN_R1="${NONHUMAN_R1}"
    KRAKEN_R2="${NONHUMAN_R2}"
  else
    log "Host removal skipped: ${SAMPLE_ID}"
    KRAKEN_R1="${CLEAN_R1}"
    KRAKEN_R2="${CLEAN_R2}"
  fi

  ############################
  # Kraken2 PE
  ############################

  if [[ -s "${KRAKEN_DONE}" && ! -s "${KREPORT}" ]]; then
    log "Kraken2 .done flag is inconsistent with report; rebuilding ${SAMPLE_ID}"
    rm -f "${KRAKEN_DONE}" "${KREPORT}"
  fi

  if [[ ! -s "${KRAKEN_DONE}" ]]; then
    stage_start "${SAMPLE_ID}" "kraken"
    log "Kraken2 PE: ${SAMPLE_ID}"

    rm -f "${KREPORT_PART}"
    if [[ "${KEEP_KRAKEN_OUTPUT}" == "1" ]]; then
      rm -f "${KOUTPUT_PART}"
    fi

    kraken2 \
      --db "${KRAKEN_DB}" \
      --threads "${KRAKEN_THREADS}" \
      --paired \
      --gzip-compressed \
      --use-names \
      --report "${KREPORT_PART}" \
      --output "${KOUTPUT_PART}" \
      "${KRAKEN_R1}" "${KRAKEN_R2}" \
      > "${LOG_DIR}/${SAMPLE_ID}.kraken2.log" 2>&1

    [[ -s "${KREPORT_PART}" ]] || \
      fail "Kraken2 report is empty for ${SAMPLE_ID}"

    mv -f "${KREPORT_PART}" "${KREPORT}"

    if [[ "${KEEP_KRAKEN_OUTPUT}" == "1" ]]; then
      [[ -s "${KOUTPUT_PART}" ]] || \
        fail "Kraken2 per-read output is empty for ${SAMPLE_ID}"
      mv -f "${KOUTPUT_PART}" "${KOUTPUT}"
    fi

    stage_finish "${SAMPLE_ID}" "kraken"
  else
    log "Completed Kraken2 report detected: ${SAMPLE_ID}"
  fi

  [[ -s "${KREPORT}" ]] || fail "Kraken2 report is empty: ${KREPORT}"

  ############################
  # Bracken
  ############################

  stage_start "${SAMPLE_ID}" "bracken"

  local LVL
  local BRACKEN_OUT
  local BRACKEN_PART
  local BRACKEN_LEVEL_DONE

  for LVL in "${BRACKEN_LEVELS[@]}"; do
    BRACKEN_OUT="${BRACKEN_DIR}/${LVL}/${SAMPLE_ID}.bracken.tsv"
    BRACKEN_PART="${BRACKEN_DIR}/${LVL}/.${SAMPLE_ID}.part.bracken.tsv"
    BRACKEN_LEVEL_DONE="${STATUS_DIR}/${SAMPLE_ID}.bracken.${LVL}.done"

    if [[ -s "${BRACKEN_LEVEL_DONE}" && ! -s "${BRACKEN_OUT}" ]]; then
      rm -f "${BRACKEN_LEVEL_DONE}"
    fi

    if [[ ! -s "${BRACKEN_LEVEL_DONE}" ]]; then
      log "Bracken level=${LVL}: ${SAMPLE_ID}"
      rm -f "${BRACKEN_PART}"

      bracken \
        -d "${KRAKEN_DB}" \
        -i "${KREPORT}" \
        -o "${BRACKEN_PART}" \
        -r "${READ_LEN}" \
        -l "${LVL}" \
        -t "${BRACKEN_THRESHOLD}" \
        > "${LOG_DIR}/${SAMPLE_ID}.bracken.${LVL}.log" 2>&1

      [[ -s "${BRACKEN_PART}" ]] || \
        fail "Bracken output is empty: ${SAMPLE_ID}, level=${LVL}"

      mv -f "${BRACKEN_PART}" "${BRACKEN_OUT}"
      mark_time "${BRACKEN_LEVEL_DONE}"
    else
      log "Completed Bracken output level=${LVL}: ${SAMPLE_ID}"
    fi
  done

  stage_finish "${SAMPLE_ID}" "bracken"

  sample_outputs_complete "${SAMPLE_ID}" || \
    fail "Final Kraken2/Bracken outputs are incomplete for ${SAMPLE_ID}"

  mark_time "${SAMPLE_DONE}"
  rm -f "${SAMPLE_RUNNING}" "${SAMPLE_FAILED}"
  ACTIVE_STAGE=""

  ############################
  # Cleanup after sample done
  ############################

  if [[ "${CLEAN_INTERMEDIATE}" == "1" ]]; then
    log "Cleaning intermediate files after successful sample completion: ${SAMPLE_ID}"
    rm -f \
      "${CLEAN_R1}" "${CLEAN_R2}" \
      "${NONHUMAN_R1}" "${NONHUMAN_R2}" \
      "${HOST_BAM}"
    rm -rf "${WORK_DIR}"
  fi

  log "Finished sample: ${SAMPLE_ID}"
}

###############################################################################
# Run samples serially
###############################################################################

log "Starting per-sample analysis..."

TOTAL_SAMPLES=$(awk -F'\t' 'NR > 1 {count++} END {print count + 0}' "${MANIFEST}")
CURRENT_SAMPLE=0

while IFS=$'\t' read -r SAMPLE_ID LAYOUT R1 R2 STATUS; do
  [[ "${SAMPLE_ID}" == "sample_id" ]] && continue
  [[ -n "${SAMPLE_ID}" ]] || continue

  CURRENT_SAMPLE=$((CURRENT_SAMPLE + 1))
  log "Progress: ${CURRENT_SAMPLE}/${TOTAL_SAMPLES} (${SAMPLE_ID})"

  ACTIVE_SAMPLE="${SAMPLE_ID}"
  process_sample "${SAMPLE_ID}" "${R1}" "${R2}"
  ACTIVE_SAMPLE=""
  ACTIVE_STAGE=""
done < "${MANIFEST}"

if [[ "${AC_STOP_AFTER_QC:-0}" == "1" ]]; then
  log "All available FASTQ pairs passed the recovered fastp preprocessing stage."
  exit 0
fi

###############################################################################
# Build completed manifest
###############################################################################

COMPLETED_MANIFEST="${MANIFEST_DIR}/manifest_completed.tsv"

python3 - <<PY
from pathlib import Path
import pandas as pd

manifest = pd.read_csv(r"${MANIFEST}", sep="\t")
status_dir = Path(r"${STATUS_DIR}")

completed = manifest[
    manifest["sample_id"].astype(str).map(
        lambda sample: (status_dir / f"{sample}.done").exists()
    )
].copy()

completed.to_csv(r"${COMPLETED_MANIFEST}", sep="\t", index=False)
print(f"[INFO] Completed samples: {len(completed)}/{len(manifest)}")

if len(completed) == 0:
    raise SystemExit("[ERROR] No completed samples are available for merging.")
PY

###############################################################################
# Merge Bracken outputs
###############################################################################

MERGE_SCRIPT="${SCRIPT_DIR}/merge_bracken_to_16s_like_tables.py"

cat > "${MERGE_SCRIPT}" <<'PYMERGE'
#!/usr/bin/env python3

from pathlib import Path
import argparse
import re
import pandas as pd

RANK_NAME = {
    "S": "species",
    "G": "genus",
    "F": "family",
    "P": "phylum",
}

RANK_PREFIX = {
    "D": "d__",
    "K": "k__",
    "P": "p__",
    "C": "c__",
    "O": "o__",
    "F": "f__",
    "G": "g__",
    "S": "s__",
}

MAIN_RANKS = ["D", "K", "P", "C", "O", "F", "G", "S"]


def clean_name(value):
    value = str(value).strip()
    return re.sub(r"\s+", "_", value)


def parse_kraken_report(path):
    lineage_map = {}
    stack = []

    with open(path, "r", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 6:
                continue

            rank_code = parts[3].strip()
            taxid = parts[4].strip()
            raw_name = parts[5]
            name = raw_name.strip()

            indent = len(raw_name) - len(raw_name.lstrip(" "))
            depth = indent // 2
            rank_simple = rank_code[0] if rank_code else ""

            while stack and stack[-1][0] >= depth:
                stack.pop()

            stack.append((depth, rank_simple, taxid, name))

            seen = {}
            for _, rank, _, tax_name in stack:
                if rank in RANK_PREFIX:
                    seen[rank] = clean_name(tax_name)

            lineage = [
                RANK_PREFIX[rank] + seen[rank]
                for rank in MAIN_RANKS
                if rank in seen
            ]

            if lineage:
                lineage_map[taxid] = ";".join(lineage)

    return lineage_map


def read_bracken(path):
    try:
        frame = pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

    required = [
        "name",
        "taxonomy_id",
        "taxonomy_lvl",
        "kraken_assigned_reads",
        "added_reads",
        "new_est_reads",
        "fraction_total_reads",
    ]

    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")

    frame["taxonomy_id"] = frame["taxonomy_id"].astype(str)
    frame["Feature ID"] = "taxid_" + frame["taxonomy_id"]
    frame["new_est_reads"] = pd.to_numeric(
        frame["new_est_reads"], errors="coerce"
    ).fillna(0)
    frame["fraction_total_reads"] = pd.to_numeric(
        frame["fraction_total_reads"], errors="coerce"
    ).fillna(0)
    return frame


def merge_tables(tables):
    merged = tables[0]
    for table in tables[1:]:
        merged = merged.merge(table, on="Feature ID", how="outer")
    return merged.fillna(0)


def write_level(
    level,
    sample_ids,
    bracken_dir,
    kraken_report_dir,
    export_dir,
    tax_dir,
):
    rank_name = RANK_NAME[level]
    level_dir = bracken_dir / level

    count_tables = []
    relative_tables = []
    taxonomy_records = {}
    long_tables = []
    lineage_map = {}
    used_samples = []

    for sample in sample_ids:
        report = kraken_report_dir / f"{sample}.kreport"
        if report.exists() and report.stat().st_size > 0:
            lineage_map.update(parse_kraken_report(report))

    for sample in sample_ids:
        bracken_file = level_dir / f"{sample}.bracken.tsv"

        if not bracken_file.exists() or bracken_file.stat().st_size == 0:
            print(f"[WARNING] Missing Bracken file: {bracken_file}")
            continue

        frame = read_bracken(bracken_file)
        if frame.empty:
            print(f"[WARNING] Empty Bracken table: {bracken_file}")
            continue

        used_samples.append(sample)
        frame["SampleID"] = sample

        count_tables.append(
            frame[["Feature ID", "new_est_reads"]].rename(
                columns={"new_est_reads": sample}
            )
        )
        relative_tables.append(
            frame[["Feature ID", "fraction_total_reads"]].rename(
                columns={"fraction_total_reads": sample}
            )
        )

        for _, row in frame.iterrows():
            feature_id = str(row["Feature ID"])
            taxid = str(row["taxonomy_id"])
            name = str(row["name"])
            taxon = lineage_map.get(taxid)

            if not taxon:
                prefix = {
                    "S": "s__",
                    "G": "g__",
                    "F": "f__",
                    "P": "p__",
                }[level]
                taxon = prefix + clean_name(name)

            taxonomy_records[feature_id] = {
                "Feature ID": feature_id,
                "Taxon": taxon,
                "Confidence": "",
                "TaxID": taxid,
                "Name": name,
                "Rank": level,
            }

        frame["Taxon"] = frame["taxonomy_id"].map(lineage_map)
        long_tables.append(frame)

    if not count_tables:
        raise RuntimeError(f"No non-empty Bracken tables found for level {level}")

    count_merged = merge_tables(count_tables)
    relative_merged = merge_tables(relative_tables)
    taxonomy = pd.DataFrame(
        list(taxonomy_records.values())
    ).sort_values(["Taxon", "Feature ID"])
    long_all = pd.concat(long_tables, ignore_index=True)

    feature_order = taxonomy["Feature ID"].tolist()

    count_merged = (
        count_merged.set_index("Feature ID")
        .reindex(feature_order)
        .fillna(0)
        .reset_index()
    )
    relative_merged = (
        relative_merged.set_index("Feature ID")
        .reindex(feature_order)
        .fillna(0)
        .reset_index()
    )

    if level == "S":
        count_out = export_dir / "feature-table.tsv"
        relative_out = export_dir / "feature-table-relative.tsv"
        taxonomy_out = tax_dir / "taxonomy.tsv"
        qiime_taxonomy_out = tax_dir / "taxonomy_qiime_import.tsv"
        long_out = export_dir / "long-table.tsv"
        sample_count_out = export_dir / "sample-by-feature-table.tsv"
        sample_relative_out = export_dir / "sample-by-feature-table-relative.tsv"
    else:
        count_out = export_dir / f"{rank_name}-feature-table.tsv"
        relative_out = export_dir / f"{rank_name}-feature-table-relative.tsv"
        taxonomy_out = tax_dir / f"{rank_name}-taxonomy.tsv"
        qiime_taxonomy_out = tax_dir / f"{rank_name}-taxonomy_qiime_import.tsv"
        long_out = export_dir / f"{rank_name}-long-table.tsv"
        sample_count_out = (
            export_dir / f"{rank_name}-sample-by-feature-table.tsv"
        )
        sample_relative_out = (
            export_dir / f"{rank_name}-sample-by-feature-table-relative.tsv"
        )

    count_merged.rename(
        columns={"Feature ID": "#OTU ID"}
    ).to_csv(count_out, sep="\t", index=False)

    relative_merged.rename(
        columns={"Feature ID": "#OTU ID"}
    ).to_csv(relative_out, sep="\t", index=False)

    taxonomy[
        ["Feature ID", "Taxon", "Confidence"]
    ].to_csv(taxonomy_out, sep="\t", index=False)

    taxonomy[
        ["Feature ID", "Taxon"]
    ].to_csv(
        qiime_taxonomy_out,
        sep="\t",
        index=False,
        header=False,
    )

    long_all.to_csv(long_out, sep="\t", index=False)

    count_sample_by_feature = count_merged.set_index("Feature ID").T
    count_sample_by_feature.index.name = "SampleID"
    count_sample_by_feature.to_csv(sample_count_out, sep="\t")

    relative_sample_by_feature = relative_merged.set_index("Feature ID").T
    relative_sample_by_feature.index.name = "SampleID"
    relative_sample_by_feature.to_csv(sample_relative_out, sep="\t")

    print(
        f"[INFO] {rank_name}: "
        f"{len(feature_order)} features across {len(used_samples)} samples"
    )
    print(f"[INFO] Wrote: {count_out}")
    print(f"[INFO] Wrote: {relative_out}")
    print(f"[INFO] Wrote: {taxonomy_out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--bracken_dir", required=True)
    parser.add_argument("--kraken_report_dir", required=True)
    parser.add_argument("--export_dir", required=True)
    parser.add_argument("--tax_dir", required=True)
    parser.add_argument(
        "--levels",
        nargs="+",
        default=["S", "G", "F", "P"],
    )
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")
    sample_ids = manifest["sample_id"].astype(str).tolist()

    bracken_dir = Path(args.bracken_dir)
    kraken_report_dir = Path(args.kraken_report_dir)
    export_dir = Path(args.export_dir)
    tax_dir = Path(args.tax_dir)

    export_dir.mkdir(parents=True, exist_ok=True)
    tax_dir.mkdir(parents=True, exist_ok=True)

    for level in args.levels:
        write_level(
            level,
            sample_ids,
            bracken_dir,
            kraken_report_dir,
            export_dir,
            tax_dir,
        )


if __name__ == "__main__":
    main()
PYMERGE

chmod +x "${MERGE_SCRIPT}"

log "Merging Bracken outputs..."

python3 "${MERGE_SCRIPT}" \
  --manifest "${COMPLETED_MANIFEST}" \
  --bracken_dir "${BRACKEN_DIR}" \
  --kraken_report_dir "${KRAKEN_REPORT_DIR}" \
  --export_dir "${EXPORT_DIR}" \
  --tax_dir "${TAX_DIR}" \
  --levels "${BRACKEN_LEVELS[@]}"

###############################################################################
# Optional BIOM and QIIME2 imports
###############################################################################

if command -v biom >/dev/null 2>&1; then
  log "Creating BIOM tables..."

  if [[ -s "${EXPORT_DIR}/feature-table.tsv" ]]; then
    biom convert \
      -i "${EXPORT_DIR}/feature-table.tsv" \
      -o "${EXPORT_DIR}/feature-table.biom" \
      --to-hdf5 \
      --table-type "OTU table" \
      || log "WARNING: BIOM conversion failed for species table."
  fi

  if [[ -s "${EXPORT_DIR}/genus-feature-table.tsv" ]]; then
    biom convert \
      -i "${EXPORT_DIR}/genus-feature-table.tsv" \
      -o "${EXPORT_DIR}/genus-feature-table.biom" \
      --to-hdf5 \
      --table-type "OTU table" \
      || log "WARNING: BIOM conversion failed for genus table."
  fi
else
  log "WARNING: biom command not found; BIOM conversion skipped."
fi

if command -v qiime >/dev/null 2>&1 && [[ -s "${EXPORT_DIR}/feature-table.biom" ]]; then
  log "Importing species table and taxonomy into QIIME2..."

  qiime tools import \
    --input-path "${EXPORT_DIR}/feature-table.biom" \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path "${EXPORT_DIR}/table.qza" \
    || log "WARNING: QIIME2 table import failed."

  if [[ -s "${EXPORT_DIR}/table.qza" ]]; then
    qiime feature-table summarize \
      --i-table "${EXPORT_DIR}/table.qza" \
      --o-visualization "${EXPORT_DIR}/table.qzv" \
      || log "WARNING: QIIME2 table summarization failed."
  fi

  qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path "${TAX_DIR}/taxonomy_qiime_import.tsv" \
    --output-path "${TAX_DIR}/taxonomy.qza" \
    || log "WARNING: QIIME2 taxonomy import failed."
else
  log "QIIME2 not found or species BIOM missing; QZA/QZV generation skipped."
fi

###############################################################################
# Final summary
###############################################################################

{
  printf 'run_end\t%s\n' "$(date -Iseconds)"
  printf 'completed_manifest\t%s\n' "${COMPLETED_MANIFEST}"
} >> "${PROJECT_INFO}"

log "=================================================================="
log "Metagenome taxonomy profiling finished."
log "Species count table: ${EXPORT_DIR}/feature-table.tsv"
log "Species relative table: ${EXPORT_DIR}/feature-table-relative.tsv"
log "Species taxonomy: ${TAX_DIR}/taxonomy.tsv"
log "Genus count table: ${EXPORT_DIR}/genus-feature-table.tsv"
log "Genus relative table: ${EXPORT_DIR}/genus-feature-table-relative.tsv"
log "Genus taxonomy: ${TAX_DIR}/genus-taxonomy.tsv"
log "Manifest: ${MANIFEST}"
log "Completed manifest: ${COMPLETED_MANIFEST}"
log "Read-length QC: ${LENGTH_QC}"
log "=================================================================="
