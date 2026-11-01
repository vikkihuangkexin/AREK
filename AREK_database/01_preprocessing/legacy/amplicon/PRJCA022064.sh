#!/usr/bin/env bash
set -euo pipefail
# Portable recovered source. Canonical resumable route: run_pipeline.py.
: "${RAW_DATA_ROOT:?set raw data root}" "${RESULT_ROOT:?set isolated output root}" "${DATABASE_ROOT:?set reference root}"
TMP_ROOT="${TMP_ROOT:-${RESULT_ROOT}/tmp}"
if [[ "${AREK_LEGACY_SOURCE_MODE:-}" != "isolated_research" ]]; then echo "Recovered historical entry; use canonical run_pipeline.py. To inspect native batching/diagnostics in a NEW output directory set AREK_LEGACY_SOURCE_MODE=isolated_research." >&2; exit 2; fi


PROJECT="${PROJECT:-PRJCA022064}"
RAW_ROOT="${RAW_ROOT:-${RAW_DATA_ROOT}}"
RESULT_ROOT="${RESULT_ROOT:-${RESULT_ROOT}}"
CLASSIFIER="${CLASSIFIER:-${DATABASE_ROOT}/silva/silva-138-99-nb-classifier.qza}"
THREADS="${THREADS:-16}"
PERC_IDENTITY="${PERC_IDENTITY:-0.97}"
MIN_MEDIAN_LEN="${MIN_MEDIAN_LEN:-100}"
QC_READS_PER_FILE="${QC_READS_PER_FILE:-3000}"
FORCE="${FORCE:-0}"

# If paired merging fails due to insufficient overlap, set FORCE_SINGLE_BEST_MATE=1 to use the longer mate as single-end input for all paired samples.
FORCE_SINGLE_BEST_MATE="${FORCE_SINGLE_BEST_MATE:-0}"

PROJECT_DIR="${RAW_ROOT}/${PROJECT}"
FASTQ_DIR="${PROJECT_DIR}/fastq"
RESULT_DIR="${RESULT_ROOT}/${PROJECT}"

MANIFEST_DIR="${RESULT_DIR}/00_manifest"
IMPORT_DIR="${RESULT_DIR}/01_import"
DEMUX_DIR="${RESULT_DIR}/02_demux"
VSEARCH_DIR="${RESULT_DIR}/03_vsearch_otu97"
EXPORT_DIR="${RESULT_DIR}/04_export"
TAX_DIR="${RESULT_DIR}/05_taxonomy"
LOG_DIR="${RESULT_DIR}/logs"

mkdir -p "$MANIFEST_DIR" "$IMPORT_DIR" "$DEMUX_DIR" "$VSEARCH_DIR" "$EXPORT_DIR" "$TAX_DIR" "$LOG_DIR"

log() {
  echo "[$(date '+%F %T')] $*" | tee -a "$LOG_DIR/run_PRJCA022064_vsearch.log"
}

if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "[ERROR] FASTQ_DIR does not exist: $FASTQ_DIR" >&2
  exit 1
fi

if [[ ! -f "$CLASSIFIER" ]]; then
  echo "[ERROR] CLASSIFIER does not exist: $CLASSIFIER" >&2
  exit 1
fi

if ! qiime vsearch --help >/dev/null 2>&1; then
  echo "[ERROR] The active QIIME2 environment is missing q2-vsearch. Install the plugin before continuing." >&2
  exit 1
fi

if [[ "$FORCE" == "1" ]]; then
  log "FORCE=1: remove previous outputs except logs"
  rm -rf "$MANIFEST_DIR" "$IMPORT_DIR" "$DEMUX_DIR" "$VSEARCH_DIR" "$EXPORT_DIR" "$TAX_DIR"
  mkdir -p "$MANIFEST_DIR" "$IMPORT_DIR" "$DEMUX_DIR" "$VSEARCH_DIR" "$EXPORT_DIR" "$TAX_DIR"
fi

cat <<EOF | tee -a "$LOG_DIR/run_PRJCA022064_vsearch.log"
============================================================
Start PRJCA022064 QIIME2 + VSEARCH workflow

PROJECT:                  $PROJECT
FASTQ_DIR:                $FASTQ_DIR
RESULT_DIR:               $RESULT_DIR
CLASSIFIER:               $CLASSIFIER
THREADS:                  $THREADS
PERC_IDENTITY:            $PERC_IDENTITY
MIN_MEDIAN_LEN:           $MIN_MEDIAN_LEN
QC_READS_PER_FILE:        $QC_READS_PER_FILE
FORCE_SINGLE_BEST_MATE:   $FORCE_SINGLE_BEST_MATE

Recognized paired suffix examples:
  *_R1.fastq.gz / *_R2.fastq.gz
  *_1.fastq.gz  / *_2.fastq.gz
  *_f1.fastq.gz / *_r2.fastq.gz
============================================================
EOF

log "Detect FASTQ layout and generate manifest"

python - <<PY
import re
import gzip
import statistics
from pathlib import Path
from collections import defaultdict, Counter

project = "$PROJECT"
fastq_dir = Path("$FASTQ_DIR")
manifest = Path("$MANIFEST_DIR/manifest.tsv")
manifest_single = Path("$MANIFEST_DIR/manifest_single.tsv")
manifest_paired = Path("$MANIFEST_DIR/manifest_paired.tsv")
layout_file = Path("$MANIFEST_DIR/read_layout.tsv")
used_file = Path("$MANIFEST_DIR/fastq_used.tsv")
all_class_file = Path("$MANIFEST_DIR/all_fastq_classification.tsv")
pair_qc_file = Path("$MANIFEST_DIR/pair_mate_length_qc.tsv")
single_qc_file = Path("$MANIFEST_DIR/single_length_qc.tsv")

min_median_len = int("$MIN_MEDIAN_LEN")
qc_reads = int("$QC_READS_PER_FILE")
force_single_best_mate = int("$FORCE_SINGLE_BEST_MATE")

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
        r'(.+)_R1_001$',
        r'(.+)_R1$',
        r'(.+)\.R1$',
        r'(.+)_1$',
        r'(.+)\.1$',
        r'(.+)_read1$',
        r'(.+)\.read1$',
        r'(.+)_f1$',
        r'(.+)\.f1$',
        r'(.+)_F1$',
        r'(.+)\.F1$',
        r'(.+)_forward$',
        r'(.+)\.forward$',
        r'(.+)_fwd$',
        r'(.+)\.fwd$'
    ]

    patterns_r2 = [
        r'(.+)_R2_001$',
        r'(.+)_R2$',
        r'(.+)\.R2$',
        r'(.+)_2$',
        r'(.+)\.2$',
        r'(.+)_read2$',
        r'(.+)\.read2$',
        r'(.+)_r2$',
        r'(.+)\.r2$',
        r'(.+)_R2$',
        r'(.+)\.R2$',
        r'(.+)_reverse$',
        r'(.+)\.reverse$',
        r'(.+)_rev$',
        r'(.+)\.rev$'
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

def fastq_qc(path, max_reads=3000):
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

with all_class_file.open("w") as out:
    out.write("filename\tstem\tclassified_read\tsample_id\tpath\n")

    for f in files:
        stem = strip_ext(f.name)
        read, sid_raw = classify_read(stem)
        sid = clean_sid(sid_raw)

        out.write(f"{f.name}\t{stem}\t{read}\t{sid}\t{f.resolve()}\n")

        if read in ("R1", "R2"):
            if read in pairs[sid]:
                raise SystemExit(f"Duplicate {read} for sample {sid}: {pairs[sid][read]} and {f.resolve()}")
            pairs[sid][read] = str(f.resolve())
        else:
            if sid in singles:
                raise SystemExit(f"Duplicate single file for sample {sid}: {singles[sid]} and {f.resolve()}")
            singles[sid] = str(f.resolve())

complete_pairs = []
incomplete_pairs = []

for sid, d in sorted(pairs.items()):
    if "R1" in d and "R2" in d:
        complete_pairs.append((sid, d["R1"], d["R2"]))
    else:
        incomplete_pairs.append((sid, d.get("R1", ""), d.get("R2", "")))

def choose_longer_mate(sid, r1, r2, q1, q2, suffix=""):
    r1_med = q1["median_len"] if q1["median_len"] != "NA" else 0
    r2_med = q2["median_len"] if q2["median_len"] != "NA" else 0

    if r1_med >= r2_med:
        return sid, r1, "single_use_R1" + suffix
    return sid, r2, "single_use_R2" + suffix

records_single = []
records_paired = []
both_short = []

with pair_qc_file.open("w") as out:
    out.write(
        "sample_id\tR1\tR2\t"
        "R1_median_len\tR2_median_len\tR1_min_len\tR2_min_len\tR1_max_len\tR2_max_len\t"
        "R1_qmin\tR1_qmax\tR2_qmin\tR2_qmax\tdecision\tR1_top_lengths\tR2_top_lengths\n"
    )

    for sid, r1, r2 in complete_pairs:
        q1 = fastq_qc(r1, qc_reads)
        q2 = fastq_qc(r2, qc_reads)

        r1_med = q1["median_len"] if q1["median_len"] != "NA" else 0
        r2_med = q2["median_len"] if q2["median_len"] != "NA" else 0
        r1_long = r1_med >= min_median_len
        r2_long = r2_med >= min_median_len

        if r1_long and r2_long and not force_single_best_mate:
            decision = "paired_merge"
            records_paired.append((sid, r1, r2))
        elif r1_long and r2_long and force_single_best_mate:
            sid2, path, decision = choose_longer_mate(sid, r1, r2, q1, q2, "_both_long_longer")
            records_single.append((sid2, path, decision))
        elif r1_long and not r2_long:
            decision = "single_use_R1"
            records_single.append((sid, r1, decision))
        elif (not r1_long) and r2_long:
            decision = "single_use_R2"
            records_single.append((sid, r2, decision))
        else:
            decision = "skip_both_short"
            both_short.append(sid)

        out.write(
            f"{sid}\t{r1}\t{r2}\t"
            f"{q1['median_len']}\t{q2['median_len']}\t"
            f"{q1['min_len']}\t{q2['min_len']}\t"
            f"{q1['max_len']}\t{q2['max_len']}\t"
            f"{q1['qmin']}\t{q1['qmax']}\t{q2['qmin']}\t{q2['qmax']}\t"
            f"{decision}\t{q1['top_lengths']}\t{q2['top_lengths']}\n"
        )

with single_qc_file.open("w") as out:
    out.write("sample_id\tfastq\tmedian_len\tmin_len\tmax_len\tqmin\tqmax\tdecision\ttop_lengths\n")

    for sid, path in sorted(singles.items()):
        q = fastq_qc(path, qc_reads)
        med = q["median_len"] if q["median_len"] != "NA" else 0

        if med >= min_median_len:
            decision = "single_use"
            records_single.append((sid, path, decision))
        else:
            decision = "skip_single_too_short"

        out.write(
            f"{sid}\t{path}\t{q['median_len']}\t{q['min_len']}\t{q['max_len']}\t"
            f"{q['qmin']}\t{q['qmax']}\t{decision}\t{q['top_lengths']}\n"
        )

# Determine final layout.
# If all usable samples are true paired and no single records, use paired merge.
# Otherwise unify everything as single-end selected biological reads.
notes = []
if records_paired and not records_single and not incomplete_pairs:
    analysis_layout = "paired_merge_then_vsearch"

    with manifest_paired.open("w") as out:
        out.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
        for sid, r1, r2 in records_paired:
            out.write(f"{sid}\t{r1}\t{r2}\n")

    manifest.write_text(manifest_paired.read_text())

    with used_file.open("w") as out:
        out.write("sample_id\tforward_fastq\treverse_fastq\tmode\n")
        for sid, r1, r2 in records_paired:
            out.write(f"{sid}\t{r1}\t{r2}\tpaired_merge\n")

    samples = len(records_paired)

else:
    # Convert paired-long records to selected single-end too, so all samples share one import type.
    # This branch handles index+read, mixed single+paired, incomplete pairs, and forced single mode.
    if records_paired:
        for sid, r1, r2 in records_paired:
            q1 = fastq_qc(r1, min(qc_reads, 1000))
            q2 = fastq_qc(r2, min(qc_reads, 1000))
            sid2, path, decision = choose_longer_mate(sid, r1, r2, q1, q2, "_paired_longer")
            records_single.append((sid2, path, decision))
        records_paired = []

    # Handle incomplete pairs by using available long mate
    for sid, r1, r2 in incomplete_pairs:
        for mate, path in [("R1", r1), ("R2", r2)]:
            if not path:
                continue
            q = fastq_qc(path, qc_reads)
            med = q["median_len"] if q["median_len"] != "NA" else 0
            if med >= min_median_len:
                records_single.append((sid, path, f"single_use_{mate}_incomplete_pair"))

    if records_single:
        analysis_layout = "single_selected_mate_vsearch"

        # avoid duplicate sample IDs after conversion
        seen = {}
        final_records = []
        for sid, path, mode in records_single:
            sid0 = sid
            if sid0 in seen:
                seen[sid0] += 1
                sid = f"{sid0}_{seen[sid0]}"
            else:
                seen[sid0] = 1
            final_records.append((sid, path, mode))

        with manifest_single.open("w") as out:
            out.write("sample-id\tabsolute-filepath\n")
            for sid, path, mode in final_records:
                out.write(f"{sid}\t{path}\n")

        manifest.write_text(manifest_single.read_text())

        with used_file.open("w") as out:
            out.write("sample_id\tfastq\tmode\n")
            for sid, path, mode in final_records:
                out.write(f"{sid}\t{path}\t{mode}\n")

        samples = len(final_records)

        notes.append(
            "Data were processed as selected single-end biological reads because at least one sample was single-end, index+read, incomplete, or FORCE_SINGLE_BEST_MATE=1."
        )
    else:
        analysis_layout = "skip_too_short"
        samples = 0
        notes.append("No usable reads with median length >= MIN_MEDIAN_LEN.")

# input_layout summary
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

with layout_file.open("w") as out:
    out.write("metric\tvalue\n")
    out.write(f"project\t{project}\n")
    out.write(f"input_layout\t{input_layout}\n")
    out.write(f"analysis_layout\t{analysis_layout}\n")
    out.write(f"samples\t{samples}\n")
    out.write(f"fastq_files_found\t{len(files)}\n")
    out.write(f"complete_pairs\t{len(complete_pairs)}\n")
    out.write(f"incomplete_pairs\t{len(incomplete_pairs)}\n")
    out.write(f"single_files\t{len(singles)}\n")
    out.write(f"skipped_both_short_pairs\t{len(both_short)}\n")
    out.write(f"min_median_len\t{min_median_len}\n")
    out.write(f"qc_reads_per_file\t{qc_reads}\n")
    out.write(f"force_single_best_mate\t{force_single_best_mate}\n")
    out.write(f"search_dir\t{fastq_dir}\n")
    if notes:
        out.write("notes\t" + " ".join(notes).replace("\t", " ") + "\n")

print(layout_file.read_text())
PY

cat "$MANIFEST_DIR/read_layout.tsv" | tee -a "$LOG_DIR/run_PRJCA022064_vsearch.log"

ANALYSIS_LAYOUT="$(awk -F'\t' '$1=="analysis_layout"{print $2}' "$MANIFEST_DIR/read_layout.tsv")"

if [[ "$ANALYSIS_LAYOUT" == "skip_too_short" ]]; then
  echo "[SKIP] $PROJECT has no usable reads." | tee "$LOG_DIR/skipped.too_short_reads.log"
  exit 0
fi

if [[ "$ANALYSIS_LAYOUT" == "paired_merge_then_vsearch" ]]; then
  RAW_DEMUX_QZA="$IMPORT_DIR/paired-end-demux.qza"
  MERGED_DEMUX_QZA="$IMPORT_DIR/merged-demux.qza"
  UNMERGED_DEMUX_QZA="$IMPORT_DIR/unmerged-demux.qza"

  log "QIIME2 import as paired-end"
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$MANIFEST_DIR/manifest.tsv" \
    --output-path "$RAW_DEMUX_QZA" \
    --input-format PairedEndFastqManifestPhred33V2

  log "demux summarize: raw paired-end"
  qiime demux summarize \
    --i-data "$RAW_DEMUX_QZA" \
    --o-visualization "$DEMUX_DIR/demux-raw-summary.qzv"

  log "VSEARCH merge-pairs"
  if qiime vsearch merge-pairs --help >/dev/null 2>&1; then
    qiime vsearch merge-pairs \
      --i-demultiplexed-seqs "$RAW_DEMUX_QZA" \
      --p-threads "$THREADS" \
      --o-merged-sequences "$MERGED_DEMUX_QZA" \
      --o-unmerged-sequences "$UNMERGED_DEMUX_QZA"
  else
    qiime vsearch join-pairs \
      --i-demultiplexed-seqs "$RAW_DEMUX_QZA" \
      --p-threads "$THREADS" \
      --o-joined-sequences "$MERGED_DEMUX_QZA"
  fi

  log "demux summarize: merged paired-end"
  qiime demux summarize \
    --i-data "$MERGED_DEMUX_QZA" \
    --o-visualization "$DEMUX_DIR/demux-merged-summary.qzv"

  DEMUX_FOR_VSEARCH="$MERGED_DEMUX_QZA"

else
  RAW_DEMUX_QZA="$IMPORT_DIR/single-end-demux.qza"

  log "QIIME2 import as single-end / selected biological mate"
  qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "$MANIFEST_DIR/manifest.tsv" \
    --output-path "$RAW_DEMUX_QZA" \
    --input-format SingleEndFastqManifestPhred33V2

  log "demux summarize: single-end / selected biological mate"
  qiime demux summarize \
    --i-data "$RAW_DEMUX_QZA" \
    --o-visualization "$DEMUX_DIR/demux-raw-summary.qzv"

  DEMUX_FOR_VSEARCH="$RAW_DEMUX_QZA"
fi

DEREP_TABLE="$VSEARCH_DIR/derep-table.qza"
DEREP_SEQS="$VSEARCH_DIR/derep-seqs.qza"
CLUSTER_TABLE="$VSEARCH_DIR/table.clustered.qza"
CLUSTER_SEQS="$VSEARCH_DIR/rep-seqs.clustered.qza"
CHIMERAS="$VSEARCH_DIR/chimeras.qza"
NONCHIMERAS="$VSEARCH_DIR/nonchimeras.qza"
UCHIME_STATS="$VSEARCH_DIR/uchime-stats.qza"
TABLE_QZA="$VSEARCH_DIR/table.qza"
REP_SEQS_QZA="$VSEARCH_DIR/rep-seqs.qza"

log "VSEARCH dereplicate sequences"
qiime vsearch dereplicate-sequences \
  --i-sequences "$DEMUX_FOR_VSEARCH" \
  --o-dereplicated-table "$DEREP_TABLE" \
  --o-dereplicated-sequences "$DEREP_SEQS"

log "VSEARCH de novo clustering at identity=$PERC_IDENTITY"
qiime vsearch cluster-features-de-novo \
  --i-table "$DEREP_TABLE" \
  --i-sequences "$DEREP_SEQS" \
  --p-perc-identity "$PERC_IDENTITY" \
  --p-threads "$THREADS" \
  --o-clustered-table "$CLUSTER_TABLE" \
  --o-clustered-sequences "$CLUSTER_SEQS"

log "VSEARCH uchime-denovo chimera detection"
qiime vsearch uchime-denovo \
  --i-table "$CLUSTER_TABLE" \
  --i-sequences "$CLUSTER_SEQS" \
  --o-chimeras "$CHIMERAS" \
  --o-nonchimeras "$NONCHIMERAS" \
  --o-stats "$UCHIME_STATS"

log "Filter chimeric features from table"
qiime feature-table filter-features \
  --i-table "$CLUSTER_TABLE" \
  --m-metadata-file "$NONCHIMERAS" \
  --o-filtered-table "$TABLE_QZA"

log "Filter chimeric features from representative sequences"
qiime feature-table filter-seqs \
  --i-data "$CLUSTER_SEQS" \
  --m-metadata-file "$NONCHIMERAS" \
  --o-filtered-data "$REP_SEQS_QZA"

log "Generate feature-table summary"
qiime feature-table summarize \
  --i-table "$TABLE_QZA" \
  --o-visualization "$VSEARCH_DIR/table.qzv" || true

log "Generate rep-seqs summary"
qiime feature-table tabulate-seqs \
  --i-data "$REP_SEQS_QZA" \
  --o-visualization "$VSEARCH_DIR/rep-seqs.qzv" || true

log "Export feature table"
rm -rf "$EXPORT_DIR"
mkdir -p "$EXPORT_DIR"

qiime tools export \
  --input-path "$TABLE_QZA" \
  --output-path "$EXPORT_DIR"

biom convert \
  -i "$EXPORT_DIR/feature-table.biom" \
  -o "$EXPORT_DIR/feature-table.tsv" \
  --to-tsv

log "Taxonomic classification"
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$REP_SEQS_QZA" \
  --p-n-jobs "$THREADS" \
  --o-classification "$TAX_DIR/taxonomy.qza"

qiime tools export \
  --input-path "$TAX_DIR/taxonomy.qza" \
  --output-path "$TAX_DIR"

log "Done: $PROJECT"

cat <<EOF
============================================================
Finished $PROJECT

Main outputs:
  $EXPORT_DIR/feature-table.tsv
  $TAX_DIR/taxonomy.tsv
  $VSEARCH_DIR/table.qza
  $VSEARCH_DIR/rep-seqs.qza

QC files:
  $MANIFEST_DIR/read_layout.tsv
  $MANIFEST_DIR/fastq_used.tsv
  $MANIFEST_DIR/all_fastq_classification.tsv
  $MANIFEST_DIR/pair_mate_length_qc.tsv
  $MANIFEST_DIR/single_length_qc.tsv
============================================================
EOF
