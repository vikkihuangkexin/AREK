#!/usr/bin/env bash
set -euo pipefail
# Portable recovered source. Canonical resumable route: run_pipeline.py.
: "${RAW_DATA_ROOT:?set raw data root}" "${RESULT_ROOT:?set isolated output root}" "${DATABASE_ROOT:?set reference root}"
TMP_ROOT="${TMP_ROOT:-${RESULT_ROOT}/tmp}"
if [[ "${AREK_LEGACY_SOURCE_MODE:-}" != "isolated_research" ]]; then echo "Recovered historical entry; use canonical run_pipeline.py. To inspect native batching/diagnostics in a NEW output directory set AREK_LEGACY_SOURCE_MODE=isolated_research." >&2; exit 2; fi


PROJECT="${1:-${PROJECT:-}}"
if [[ -z "$PROJECT" ]]; then
  echo "Usage: PROJECT=PRJNAxxxx bash $0" >&2
  echo "   or: bash $0 PRJNAxxxx" >&2
  exit 1
fi

RAW_ROOT="${RAW_ROOT:-${RAW_DATA_ROOT}}"
RESULT_ROOT="${RESULT_ROOT:-${RESULT_ROOT}}"
CLASSIFIER="${CLASSIFIER:-${DATABASE_ROOT}/silva/silva-138-99-nb-classifier.qza}"

THREADS="${THREADS:-16}"
BATCH_SIZE="${BATCH_SIZE:-50}"
PERC_IDENTITY="${PERC_IDENTITY:-0.97}"
MIN_MEDIAN_LEN="${MIN_MEDIAN_LEN:-100}"
QC_READS_PER_FILE="${QC_READS_PER_FILE:-1000}"
BOTH_LONG_MATE="${BOTH_LONG_MATE:-R1}"
FORCE="${FORCE:-0}"

PROJECT_DIR="${RAW_ROOT}/${PROJECT}"
FASTQ_DIR="${PROJECT_DIR}/fastq"
RESULT_DIR="${RESULT_ROOT}/${PROJECT}"

MANIFEST_DIR="${RESULT_DIR}/00_manifest"
BATCH_ROOT="${RESULT_DIR}/01_batches"
MERGE_DIR="${RESULT_DIR}/02_merged_derep"
VSEARCH_DIR="${RESULT_DIR}/03_vsearch_otu97"
EXPORT_DIR="${RESULT_DIR}/04_export"
TAX_DIR="${RESULT_DIR}/05_taxonomy"
LOG_DIR="${RESULT_DIR}/logs"

TMP_BASE="${TMP_BASE:-${TMP_ROOT}/${PROJECT}}"
mkdir -p "$TMP_BASE"
export TMPDIR="$TMP_BASE"
export TMP="$TMP_BASE"
export TEMP="$TMP_BASE"

mkdir -p "$MANIFEST_DIR" "$BATCH_ROOT" "$MERGE_DIR" "$VSEARCH_DIR" "$EXPORT_DIR" "$TAX_DIR" "$LOG_DIR"

LOG_FILE="$LOG_DIR/run_${PROJECT}_batch_vsearch.log"

log() {
  echo "[$(date '+%F %T')] $*" | tee -a "$LOG_FILE"
}

if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "[ERROR] FASTQ_DIR does not exist: $FASTQ_DIR" >&2
  exit 1
fi

if [[ ! -f "$CLASSIFIER" ]]; then
  echo "[ERROR] CLASSIFIER does not exist: $CLASSIFIER" >&2
  exit 1
fi

if ! qiime --help | grep -q "vsearch"; then
  echo "[ERROR] 当前 QIIME2 环境没有 q2-vsearch 插件。" >&2
  exit 1
fi

if [[ "$FORCE" == "1" ]]; then
  log "FORCE=1: remove previous outputs except logs"
  rm -rf "$MANIFEST_DIR" "$BATCH_ROOT" "$MERGE_DIR" "$VSEARCH_DIR" "$EXPORT_DIR" "$TAX_DIR"
  mkdir -p "$MANIFEST_DIR" "$BATCH_ROOT" "$MERGE_DIR" "$VSEARCH_DIR" "$EXPORT_DIR" "$TAX_DIR"
fi

cat <<EOF | tee -a "$LOG_FILE"
============================================================
Start batch QIIME2 + VSEARCH workflow

PROJECT:             $PROJECT
FASTQ_DIR:           $FASTQ_DIR
RESULT_DIR:          $RESULT_DIR
CLASSIFIER:          $CLASSIFIER
THREADS:             $THREADS
BATCH_SIZE:          $BATCH_SIZE
PERC_IDENTITY:       $PERC_IDENTITY
MIN_MEDIAN_LEN:      $MIN_MEDIAN_LEN
QC_READS_PER_FILE:   $QC_READS_PER_FILE
BOTH_LONG_MATE:      $BOTH_LONG_MATE
TMPDIR:              $TMPDIR

Strategy:
  1. Detect FASTQ layout.
  2. For each sample, select one biological read.
  3. Split selected reads into batches.
  4. Per batch: import + dereplicate only.
  5. Merge all dereplicated tables/sequences.
  6. Global 97% OTU clustering.
  7. Chimera removal, taxonomy, final feature table.
============================================================
EOF

MANIFEST_ALL="$MANIFEST_DIR/manifest_all_selected_single.tsv"

if [[ ! -f "$MANIFEST_ALL" ]]; then
  log "Step 1: detect FASTQ layout and build selected single-end manifest"

  python - <<PY
import re
import gzip
import statistics
from pathlib import Path
from collections import defaultdict, Counter

project = "$PROJECT"
fastq_dir = Path("$FASTQ_DIR")
manifest_all = Path("$MANIFEST_ALL")
fastq_used = Path("$MANIFEST_DIR/fastq_used.tsv")
all_class = Path("$MANIFEST_DIR/all_fastq_classification.tsv")
pair_qc = Path("$MANIFEST_DIR/pair_mate_length_qc.tsv")
single_qc = Path("$MANIFEST_DIR/single_length_qc.tsv")
layout = Path("$MANIFEST_DIR/read_layout.tsv")

min_median_len = int("$MIN_MEDIAN_LEN")
qc_reads = int("$QC_READS_PER_FILE")
both_long_mate = "$BOTH_LONG_MATE"

files = sorted(
    list(fastq_dir.glob("*.fastq.gz")) +
    list(fastq_dir.glob("*.fq.gz")) +
    list(fastq_dir.glob("*.fastq")) +
    list(fastq_dir.glob("*.fq"))
)

if not files:
    raise SystemExit(f"No FASTQ files found in {fastq_dir}")

def strip_ext(name):
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(ext):
            return name[:-len(ext)]
    return name

def clean_sid(s):
    s = re.sub(r'[^A-Za-z0-9_.-]+', '_', s)
    s = re.sub(r'_+', '_', s).strip("_")
    return s or "sample"

def classify_read(stem):
    patterns_r1 = [
        r'(.+)_R1_001$', r'(.+)_R1$', r'(.+)\\.R1$',
        r'(.+)_1$', r'(.+)\\.1$',
        r'(.+)_read1$', r'(.+)\\.read1$',
        r'(.+)_f1$', r'(.+)\\.f1$',
        r'(.+)_forward$', r'(.+)\\.forward$',
        r'(.+)_fwd$', r'(.+)\\.fwd$'
    ]
    patterns_r2 = [
        r'(.+)_R2_001$', r'(.+)_R2$', r'(.+)\\.R2$',
        r'(.+)_2$', r'(.+)\\.2$',
        r'(.+)_read2$', r'(.+)\\.read2$',
        r'(.+)_r2$', r'(.+)\\.r2$',
        r'(.+)_reverse$', r'(.+)\\.reverse$',
        r'(.+)_rev$', r'(.+)\\.rev$'
    ]

    for p in patterns_r1:
        m = re.match(p, stem, flags=re.IGNORECASE)
        if m:
            return "R1", m.group(1)

    for p in patterns_r2:
        m = re.match(p, stem, flags=re.IGNORECASE)
        if m:
            return "R2", m.group(1)

    return "single", stem

def fastq_qc(path, max_reads=1000):
    opener = gzip.open if str(path).endswith(".gz") else open
    lens = []
    qmins = []
    qmaxs = []
    n = 0

    with opener(path, "rt", errors="replace") as handle:
        while True:
            h = handle.readline()
            if not h:
                break
            seq = handle.readline().strip()
            handle.readline()
            qual = handle.readline().strip()
            if not qual:
                break

            lens.append(len(seq))
            if qual:
                qs = [ord(c) for c in qual]
                qmins.append(min(qs))
                qmaxs.append(max(qs))

            n += 1
            if n >= max_reads:
                break

    if not lens:
        return {
            "n_reads": 0,
            "min_len": "NA",
            "median_len": "NA",
            "max_len": "NA",
            "qmin": "NA",
            "qmax": "NA",
            "top_lengths": ""
        }

    return {
        "n_reads": len(lens),
        "min_len": min(lens),
        "median_len": float(statistics.median(lens)),
        "max_len": max(lens),
        "qmin": min(qmins) if qmins else "NA",
        "qmax": max(qmaxs) if qmaxs else "NA",
        "top_lengths": ";".join([f"{k}:{v}" for k, v in Counter(lens).most_common(10)])
    }

pairs = defaultdict(dict)
singles = {}

with all_class.open("w") as out:
    out.write("filename\\tstem\\tclassified_read\\tsample_id\\tpath\\n")
    for f in files:
        stem = strip_ext(f.name)
        read, sid_raw = classify_read(stem)
        sid = clean_sid(sid_raw)
        out.write(f"{f.name}\\t{stem}\\t{read}\\t{sid}\\t{f.resolve()}\\n")

        if read in ("R1", "R2"):
            if read in pairs[sid]:
                raise SystemExit(f"Duplicate {read} for sample {sid}: {pairs[sid][read]} and {f.resolve()}")
            pairs[sid][read] = str(f.resolve())
        else:
            if sid in singles:
                raise SystemExit(f"Duplicate single for sample {sid}: {singles[sid]} and {f.resolve()}")
            singles[sid] = str(f.resolve())

complete_pairs = []
incomplete_pairs = []

for sid, d in sorted(pairs.items()):
    if "R1" in d and "R2" in d:
        complete_pairs.append((sid, d["R1"], d["R2"]))
    else:
        incomplete_pairs.append((sid, d.get("R1", ""), d.get("R2", "")))

records = []
skipped_single = []
skipped_pair = []

def choose_pair(sid, r1, r2, q1, q2):
    r1_med = q1["median_len"] if q1["median_len"] != "NA" else 0
    r2_med = q2["median_len"] if q2["median_len"] != "NA" else 0

    r1_long = r1_med >= min_median_len
    r2_long = r2_med >= min_median_len

    if r1_long and not r2_long:
        return sid, r1, "single_use_R1_R2_short"
    if r2_long and not r1_long:
        return sid, r2, "single_use_R2_R1_short"
    if r1_long and r2_long:
        if both_long_mate.lower() == "longer":
            if r2_med > r1_med:
                return sid, r2, "single_use_R2_both_long_longer"
            return sid, r1, "single_use_R1_both_long_longer"
        return sid, r1, "single_use_R1_both_long_default"
    return None

with pair_qc.open("w") as out:
    out.write(
        "sample_id\\tR1\\tR2\\t"
        "R1_median_len\\tR2_median_len\\tR1_min_len\\tR2_min_len\\tR1_max_len\\tR2_max_len\\t"
        "R1_qmin\\tR1_qmax\\tR2_qmin\\tR2_qmax\\tdecision\\tR1_top_lengths\\tR2_top_lengths\\n"
    )

    for sid, r1, r2 in complete_pairs:
        q1 = fastq_qc(r1, qc_reads)
        q2 = fastq_qc(r2, qc_reads)
        chosen = choose_pair(sid, r1, r2, q1, q2)

        if chosen is None:
            decision = "skip_both_short"
            skipped_pair.append(sid)
        else:
            sid2, path, decision = chosen
            records.append((sid2, path, decision))

        out.write(
            f"{sid}\\t{r1}\\t{r2}\\t"
            f"{q1['median_len']}\\t{q2['median_len']}\\t"
            f"{q1['min_len']}\\t{q2['min_len']}\\t"
            f"{q1['max_len']}\\t{q2['max_len']}\\t"
            f"{q1['qmin']}\\t{q1['qmax']}\\t{q2['qmin']}\\t{q2['qmax']}\\t"
            f"{decision}\\t{q1['top_lengths']}\\t{q2['top_lengths']}\\n"
        )

for sid, r1, r2 in incomplete_pairs:
    for mate, path in [("R1", r1), ("R2", r2)]:
        if not path:
            continue
        q = fastq_qc(path, qc_reads)
        med = q["median_len"] if q["median_len"] != "NA" else 0
        if med >= min_median_len:
            records.append((sid, path, f"single_use_{mate}_incomplete_pair"))
        else:
            skipped_pair.append(sid)

with single_qc.open("w") as out:
    out.write("sample_id\\tfastq\\tmedian_len\\tmin_len\\tmax_len\\tqmin\\tqmax\\tdecision\\ttop_lengths\\n")

    for sid, path in sorted(singles.items()):
        q = fastq_qc(path, qc_reads)
        med = q["median_len"] if q["median_len"] != "NA" else 0

        if med >= min_median_len:
            decision = "single_use"
            records.append((sid, path, decision))
        else:
            decision = "skip_single_too_short"
            skipped_single.append(sid)

        out.write(
            f"{sid}\\t{path}\\t{q['median_len']}\\t{q['min_len']}\\t{q['max_len']}\\t"
            f"{q['qmin']}\\t{q['qmax']}\\t{decision}\\t{q['top_lengths']}\\n"
        )

seen = {}
final_records = []
for sid, path, mode in records:
    base = sid
    if base in seen:
        seen[base] += 1
        sid = f"{base}_{seen[base]}"
    else:
        seen[base] = 1
    final_records.append((sid, path, mode))

if not final_records:
    raise SystemExit("No usable FASTQ reads after length QC and mate selection.")

with manifest_all.open("w") as out:
    out.write("sample-id\\tabsolute-filepath\\n")
    for sid, path, mode in final_records:
        out.write(f"{sid}\\t{path}\\n")

with fastq_used.open("w") as out:
    out.write("sample_id\\tfastq\\tmode\\n")
    for sid, path, mode in final_records:
        out.write(f"{sid}\\t{path}\\t{mode}\\n")

if complete_pairs and singles:
    input_layout = "mixed"
elif complete_pairs:
    input_layout = "paired"
elif singles:
    input_layout = "single"
elif incomplete_pairs:
    input_layout = "incomplete_pairs"
else:
    input_layout = "unknown"

with layout.open("w") as out:
    out.write("metric\\tvalue\\n")
    out.write(f"project\\t{project}\\n")
    out.write(f"input_layout\\t{input_layout}\\n")
    out.write("analysis_layout\\tbatch_selected_single_vsearch\\n")
    out.write(f"fastq_files_found\\t{len(files)}\\n")
    out.write(f"complete_pairs\\t{len(complete_pairs)}\\n")
    out.write(f"incomplete_pairs\\t{len(incomplete_pairs)}\\n")
    out.write(f"single_files\\t{len(singles)}\\n")
    out.write(f"samples_used\\t{len(final_records)}\\n")
    out.write(f"skipped_single_too_short\\t{len(skipped_single)}\\n")
    out.write(f"skipped_pair_both_short\\t{len(skipped_pair)}\\n")
    out.write(f"min_median_len\\t{min_median_len}\\n")
    out.write(f"qc_reads_per_file\\t{qc_reads}\\n")
    out.write(f"both_long_mate\\t{both_long_mate}\\n")
    out.write(f"search_dir\\t{fastq_dir}\\n")

print(layout.read_text())
PY

  cat "$MANIFEST_DIR/read_layout.tsv" | tee -a "$LOG_FILE"
else
  log "Manifest already exists, skip Step 1: $MANIFEST_ALL"
fi

log "Step 2: split manifest into batches"

BATCH_DONE="$BATCH_ROOT/batches.done"

if [[ ! -f "$BATCH_DONE" ]]; then
  rm -rf "$BATCH_ROOT"/batch_*.dir "$BATCH_ROOT"/batch_*
  mkdir -p "$BATCH_ROOT"

  N_SAMPLES="$(tail -n +2 "$MANIFEST_ALL" | wc -l)"
  N_BATCHES="$(( (N_SAMPLES + BATCH_SIZE - 1) / BATCH_SIZE ))"

  log "Samples used: $N_SAMPLES"
  log "Batch size:   $BATCH_SIZE"
  log "N batches:    $N_BATCHES"

  tail -n +2 "$MANIFEST_ALL" | split -l "$BATCH_SIZE" -d -a 4 - "$BATCH_ROOT/batch_"

  for f in "$BATCH_ROOT"/batch_*; do
    [[ -f "$f" ]] || continue
    bname="$(basename "$f")"
    bdir="$BATCH_ROOT/${bname}.dir"
    mkdir -p "$bdir"
    {
      echo -e "sample-id\tabsolute-filepath"
      cat "$f"
    } > "$bdir/manifest.tsv"
    mv "$f" "$bdir/samples.tsv"
  done

  date > "$BATCH_DONE"
else
  log "Batch split already exists, skip Step 2"
fi

log "Step 3: per-batch import + dereplicate"

mapfile -t BATCH_DIRS < <(find "$BATCH_ROOT" -maxdepth 1 -type d -name 'batch_*.dir' | sort)

if [[ "${#BATCH_DIRS[@]}" -eq 0 ]]; then
  echo "[ERROR] No batch dirs found under $BATCH_ROOT" >&2
  exit 1
fi

for bdir in "${BATCH_DIRS[@]}"; do
  bname="$(basename "$bdir" .dir)"
  log "Batch $bname: start"

  mkdir -p "$bdir/import" "$bdir/demux" "$bdir/derep"

  DEMUX_QZA="$bdir/import/single-end-demux.qza"
  DEREP_TABLE="$bdir/derep/derep-table.qza"
  DEREP_SEQS="$bdir/derep/derep-seqs.qza"
  DONE_FILE="$bdir/derep.done"

  if [[ -f "$DONE_FILE" && -f "$DEREP_TABLE" && -f "$DEREP_SEQS" ]]; then
    log "Batch $bname: already done, skip"
    continue
  fi

  if [[ ! -f "$DEMUX_QZA" ]]; then
    qiime tools import \
      --type 'SampleData[SequencesWithQuality]' \
      --input-path "$bdir/manifest.tsv" \
      --output-path "$DEMUX_QZA" \
      --input-format SingleEndFastqManifestPhred33V2
  else
    log "Batch $bname: import exists, skip"
  fi

  if [[ ! -f "$bdir/demux/demux-summary.qzv" ]]; then
    qiime demux summarize \
      --i-data "$DEMUX_QZA" \
      --o-visualization "$bdir/demux/demux-summary.qzv" || true
  fi

  qiime vsearch dereplicate-sequences \
    --i-sequences "$DEMUX_QZA" \
    --o-dereplicated-table "$DEREP_TABLE" \
    --o-dereplicated-sequences "$DEREP_SEQS"

  date > "$DONE_FILE"
  log "Batch $bname: done"
done

log "Step 4: merge all batch dereplicated tables and sequences"

MERGED_DEREP_TABLE="$MERGE_DIR/merged-derep-table.qza"
MERGED_DEREP_SEQS="$MERGE_DIR/merged-derep-seqs.qza"
MERGE_DONE="$MERGE_DIR/merge.done"

if [[ ! -f "$MERGE_DONE" || ! -f "$MERGED_DEREP_TABLE" || ! -f "$MERGED_DEREP_SEQS" ]]; then
  rm -f "$MERGED_DEREP_TABLE" "$MERGED_DEREP_SEQS" "$MERGE_DIR"/tmp-merged-*.qza

  first=1
  for bdir in "${BATCH_DIRS[@]}"; do
    table="$bdir/derep/derep-table.qza"
    seqs="$bdir/derep/derep-seqs.qza"

    if [[ ! -f "$table" || ! -f "$seqs" ]]; then
      echo "[ERROR] Missing derep output in $bdir" >&2
      exit 1
    fi

    if [[ "$first" == "1" ]]; then
      cp "$table" "$MERGED_DEREP_TABLE"
      cp "$seqs" "$MERGED_DEREP_SEQS"
      first=0
    else
      tmp_table="$MERGE_DIR/tmp-merged-table.qza"
      tmp_seqs="$MERGE_DIR/tmp-merged-seqs.qza"

      qiime feature-table merge \
        --i-tables "$MERGED_DEREP_TABLE" \
        --i-tables "$table" \
        --o-merged-table "$tmp_table"

      mv "$tmp_table" "$MERGED_DEREP_TABLE"

      qiime feature-table merge-seqs \
        --i-data "$MERGED_DEREP_SEQS" \
        --i-data "$seqs" \
        --o-merged-data "$tmp_seqs"

      mv "$tmp_seqs" "$MERGED_DEREP_SEQS"
    fi
  done

  date > "$MERGE_DONE"
else
  log "Merged derep outputs already exist, skip Step 4"
fi

log "Step 5: global 97% OTU clustering"

CLUSTER_TABLE="$VSEARCH_DIR/table.clustered.qza"
CLUSTER_SEQS="$VSEARCH_DIR/rep-seqs.clustered.qza"
CHIMERAS="$VSEARCH_DIR/chimeras.qza"
NONCHIMERAS="$VSEARCH_DIR/nonchimeras.qza"
UCHIME_STATS="$VSEARCH_DIR/uchime-stats.qza"
TABLE_QZA="$VSEARCH_DIR/table.qza"
REP_SEQS_QZA="$VSEARCH_DIR/rep-seqs.qza"

if [[ ! -f "$CLUSTER_TABLE" || ! -f "$CLUSTER_SEQS" ]]; then
  qiime vsearch cluster-features-de-novo \
    --i-table "$MERGED_DEREP_TABLE" \
    --i-sequences "$MERGED_DEREP_SEQS" \
    --p-perc-identity "$PERC_IDENTITY" \
    --p-threads "$THREADS" \
    --o-clustered-table "$CLUSTER_TABLE" \
    --o-clustered-sequences "$CLUSTER_SEQS"
else
  log "Clustered outputs exist, skip global clustering"
fi

log "Step 6: chimera detection and filtering"

if [[ ! -f "$NONCHIMERAS" ]]; then
  qiime vsearch uchime-denovo \
    --i-table "$CLUSTER_TABLE" \
    --i-sequences "$CLUSTER_SEQS" \
    --o-chimeras "$CHIMERAS" \
    --o-nonchimeras "$NONCHIMERAS" \
    --o-stats "$UCHIME_STATS"
else
  log "UCHIME outputs exist, skip"
fi

if [[ ! -f "$TABLE_QZA" ]]; then
  qiime feature-table filter-features \
    --i-table "$CLUSTER_TABLE" \
    --m-metadata-file "$NONCHIMERAS" \
    --o-filtered-table "$TABLE_QZA"
fi

if [[ ! -f "$REP_SEQS_QZA" ]]; then
  qiime feature-table filter-seqs \
    --i-data "$CLUSTER_SEQS" \
    --m-metadata-file "$NONCHIMERAS" \
    --o-filtered-data "$REP_SEQS_QZA"
fi

log "Step 7: summarize final table"

if [[ ! -f "$VSEARCH_DIR/table.qzv" ]]; then
  qiime feature-table summarize \
    --i-table "$TABLE_QZA" \
    --o-visualization "$VSEARCH_DIR/table.qzv" || true
fi

if [[ ! -f "$VSEARCH_DIR/rep-seqs.qzv" ]]; then
  qiime feature-table tabulate-seqs \
    --i-data "$REP_SEQS_QZA" \
    --o-visualization "$VSEARCH_DIR/rep-seqs.qzv" || true
fi

log "Step 8: export final feature table"

if [[ ! -f "$EXPORT_DIR/feature-table.tsv" ]]; then
  rm -rf "$EXPORT_DIR"
  mkdir -p "$EXPORT_DIR"

  qiime tools export \
    --input-path "$TABLE_QZA" \
    --output-path "$EXPORT_DIR"

  biom convert \
    -i "$EXPORT_DIR/feature-table.biom" \
    -o "$EXPORT_DIR/feature-table.tsv" \
    --to-tsv
else
  log "feature-table.tsv exists, skip export"
fi

log "Step 9: taxonomy classification"

if [[ ! -f "$TAX_DIR/taxonomy.tsv" ]]; then
  rm -rf "$TAX_DIR"
  mkdir -p "$TAX_DIR"

  qiime feature-classifier classify-sklearn \
    --i-classifier "$CLASSIFIER" \
    --i-reads "$REP_SEQS_QZA" \
    --p-n-jobs "$THREADS" \
    --o-classification "$TAX_DIR/taxonomy.qza"

  qiime tools export \
    --input-path "$TAX_DIR/taxonomy.qza" \
    --output-path "$TAX_DIR"
else
  log "taxonomy.tsv exists, skip taxonomy"
fi

log "Done: $PROJECT"

cat <<EOF
============================================================
Finished batch workflow: $PROJECT

Final outputs:
  $EXPORT_DIR/feature-table.tsv
  $TAX_DIR/taxonomy.tsv
  $VSEARCH_DIR/table.qza
  $VSEARCH_DIR/rep-seqs.qza

QC / manifest:
  $MANIFEST_DIR/read_layout.tsv
  $MANIFEST_DIR/fastq_used.tsv
  $MANIFEST_DIR/all_fastq_classification.tsv
  $MANIFEST_DIR/pair_mate_length_qc.tsv
  $MANIFEST_DIR/single_length_qc.tsv

Batches:
  $BATCH_ROOT

Merged derep:
  $MERGE_DIR/merged-derep-table.qza
  $MERGE_DIR/merged-derep-seqs.qza
============================================================
EOF
