#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "--help" ]]; then
  echo "Recovered source adapter. Read configs/source_adapters.json and INPUT_REQUIREMENTS.md."
  exit 0
fi

PROJECT="${1:-${PROJECT:-}}"
if [[ -z "$PROJECT" ]]; then
  echo "Usage: bash $0 PRJNAxxxx" >&2
  echo "Or: PROJECT=PRJNAxxxx bash $0" >&2
  exit 1
fi

RAW_ROOT="${RAW_ROOT:?Set RAW_ROOT to the original raw-read directory}"
RESULT_ROOT="${RESULT_ROOT:?Set RESULT_ROOT to an isolated output directory}"
THREADS="${THREADS:-16}"
CLASSIFIER="${CLASSIFIER:?Set CLASSIFIER to the recorded SILVA classifier}"
PERC_IDENTITY="${PERC_IDENTITY:-0.97}"
MIN_MEDIAN_LEN="${MIN_MEDIAN_LEN:-100}"
QC_READS_PER_FILE="${QC_READS_PER_FILE:-3000}"
FORCE="${FORCE:-0}"

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
  echo "[$(date '+%F %T')] $*" | tee -a "$LOG_DIR/run_qiime2_mate_aware_vsearch.log"
}

if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "[ERROR] FASTQ_DIR does not exist: $FASTQ_DIR" >&2
  exit 1
fi

if [[ ! -f "$CLASSIFIER" ]]; then
  echo "[ERROR] CLASSIFIER does not exist: $CLASSIFIER" >&2
  exit 1
fi

if [[ "$FORCE" == "1" ]]; then
  log "FORCE=1: removing previous QIIME2 outputs except logs"
  rm -rf "$MANIFEST_DIR" "$IMPORT_DIR" "$DEMUX_DIR" "$VSEARCH_DIR" "$EXPORT_DIR" "$TAX_DIR"
  mkdir -p "$MANIFEST_DIR" "$IMPORT_DIR" "$DEMUX_DIR" "$VSEARCH_DIR" "$EXPORT_DIR" "$TAX_DIR"
fi

cat <<EOF2 | tee -a "$LOG_DIR/run_qiime2_mate_aware_vsearch.log"
============================================================
[$(date '+%F %T')] Start mate-aware QIIME2 VSEARCH workflow: $PROJECT
RAW_ROOT:           $RAW_ROOT
PROJECT_DIR:        $PROJECT_DIR
FASTQ_DIR:          $FASTQ_DIR
RESULT_DIR:         $RESULT_DIR
CLASSIFIER:         $CLASSIFIER
CONDA_ENV:          ${CONDA_DEFAULT_ENV:-NA}
MODE:
  single FASTQ          -> single-end VSEARCH
  true paired FASTQ     -> merge-pairs, then VSEARCH
  index+read pairs      -> use the long biological mate as single-end
  mixed single+paired   -> convert all usable samples to single-end by selecting biological read
PERC_IDENTITY:      $PERC_IDENTITY
MIN_MEDIAN_LEN:     $MIN_MEDIAN_LEN
QC_READS_PER_FILE:  $QC_READS_PER_FILE
THREADS:            $THREADS
============================================================
EOF2

log "Detect FASTQ layout, inspect mate lengths, and generate final manifest"

python - <<PY
import os, re, gzip, statistics
from pathlib import Path
from collections import defaultdict, Counter

project = "$PROJECT"
fastq_dir = Path("$FASTQ_DIR")
manifest = Path("$MANIFEST_DIR/manifest.tsv")
manifest_single = Path("$MANIFEST_DIR/manifest_single.tsv")
manifest_paired = Path("$MANIFEST_DIR/manifest_paired.tsv")
layout_file = Path("$MANIFEST_DIR/read_layout.tsv")
used_file = Path("$MANIFEST_DIR/fastq_used.tsv")
pair_qc_file = Path("$MANIFEST_DIR/pair_mate_length_qc.tsv")
single_qc_file = Path("$MANIFEST_DIR/single_length_qc.tsv")

min_median_len = int("$MIN_MEDIAN_LEN")
qc_reads = int("$QC_READS_PER_FILE")

files = sorted(list(fastq_dir.glob("*.fastq.gz")) + list(fastq_dir.glob("*.fq.gz")) +
               list(fastq_dir.glob("*.fastq")) + list(fastq_dir.glob("*.fq")))

if not files:
    raise SystemExit(f"No FASTQ files found in {fastq_dir}")

def strip_ext(name):
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(ext):
            return name[:-len(ext)]
    return name

def classify_read(stem):
    patterns_r1 = [
        r'(.+)_R1_001$', r'(.+)_R1$', r'(.+)\.R1$',
        r'(.+)_1$', r'(.+)\.1$',
        r'(.+)_read1$', r'(.+)\.read1$'
    ]
    patterns_r2 = [
        r'(.+)_R2_001$', r'(.+)_R2$', r'(.+)\.R2$',
        r'(.+)_2$', r'(.+)\.2$',
        r'(.+)_read2$', r'(.+)\.read2$'
    ]
    for p in patterns_r1:
        m = re.match(p, stem)
        if m:
            return "R1", m.group(1)
    for p in patterns_r2:
        m = re.match(p, stem)
        if m:
            return "R2", m.group(1)
    return "single", stem

def clean_sid(s):
    s = re.sub(r'[^A-Za-z0-9_.-]+', '_', s)
    s = re.sub(r'_+', '_', s).strip("_")
    return s or "sample"

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
        return {"n_reads": 0, "min_len": "NA", "median_len": "NA", "max_len": "NA", "qmin": "NA", "qmax": "NA", "top_lengths": ""}
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

for f in files:
    stem = strip_ext(f.name)
    read, sid = classify_read(stem)
    sid = clean_sid(sid)
    if read in ("R1", "R2"):
        pairs[sid][read] = str(f.resolve())
    else:
        singles[sid] = str(f.resolve())

complete_pairs = []
incomplete_pairs = []
for sid, d in sorted(pairs.items()):
    if "R1" in d and "R2" in d:
        complete_pairs.append((sid, d["R1"], d["R2"]))
    else:
        incomplete_pairs.append((sid, d.get("R1", ""), d.get("R2", "")))

analysis_layout = None
input_layout = None
samples = 0
notes = []

if complete_pairs and not singles and not incomplete_pairs:
    input_layout = "paired"
    both_long = []
    single_best = []
    both_short = []

    with pair_qc_file.open("w") as out:
        out.write("sample_id\tR1\tR2\tR1_median_len\tR2_median_len\tR1_min_len\tR2_min_len\tR1_max_len\tR2_max_len\tR1_qmin\tR1_qmax\tR2_qmin\tR2_qmax\tdecision\n")
        for sid, r1, r2 in complete_pairs:
            q1 = fastq_qc(r1, qc_reads)
            q2 = fastq_qc(r2, qc_reads)
            r1_med = q1["median_len"] if q1["median_len"] != "NA" else 0
            r2_med = q2["median_len"] if q2["median_len"] != "NA" else 0
            r1_long = r1_med >= min_median_len
            r2_long = r2_med >= min_median_len
            if r1_long and r2_long:
                decision = "paired_merge"
                both_long.append((sid, r1, r2, r1_med, r2_med))
            elif r1_long and not r2_long:
                decision = "single_use_R1"
                single_best.append((sid, r1, decision))
            elif (not r1_long) and r2_long:
                decision = "single_use_R2"
                single_best.append((sid, r2, decision))
            else:
                decision = "skip_both_short"
                both_short.append((sid, r1, r2))
            out.write(
                f"{sid}\t{r1}\t{r2}\t{q1['median_len']}\t{q2['median_len']}\t"
                f"{q1['min_len']}\t{q2['min_len']}\t{q1['max_len']}\t{q2['max_len']}\t"
                f"{q1['qmin']}\t{q1['qmax']}\t{q2['qmin']}\t{q2['qmax']}\t{decision}\n"
            )

    if len(both_long) == len(complete_pairs):
        analysis_layout = "paired_merge_then_vsearch"
        with manifest_paired.open("w") as out:
            out.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
            for sid, r1, r2, r1_med, r2_med in both_long:
                out.write(f"{sid}\t{r1}\t{r2}\n")
        manifest.write_text(manifest_paired.read_text())
        with used_file.open("w") as out:
            out.write("sample_id\tforward_fastq\treverse_fastq\tmode\n")
            for sid, r1, r2, r1_med, r2_med in both_long:
                out.write(f"{sid}\t{r1}\t{r2}\tpaired_merge\n")
        samples = len(both_long)
    else:
        analysis_layout = "single_best_mate_vsearch"
        # For pairs where both mates are long, use the longer mate to keep one unified single-end workflow.
        for sid, r1, r2, r1_med, r2_med in both_long:
            if r2_med > r1_med:
                single_best.append((sid, r2, "single_use_R2_both_long_longer"))
            else:
                single_best.append((sid, r1, "single_use_R1_both_long_longer"))
        if not single_best:
            analysis_layout = "skip_too_short"
            notes.append("All paired samples had both mates shorter than MIN_MEDIAN_LEN.")
        else:
            with manifest_single.open("w") as out:
                out.write("sample-id\tabsolute-filepath\n")
                for sid, path, decision in single_best:
                    out.write(f"{sid}\t{path}\n")
            manifest.write_text(manifest_single.read_text())
            with used_file.open("w") as out:
                out.write("sample_id\tfastq\tmode\n")
                for sid, path, decision in single_best:
                    out.write(f"{sid}\t{path}\t{decision}\n")
            samples = len(single_best)
            notes.append(f"Not all pairs are true long-long amplicon pairs. Using long biological mate as single-end: usable={len(single_best)}, skipped_both_short={len(both_short)}.")

elif singles and not complete_pairs and not incomplete_pairs:
    input_layout = "single"
    records = []
    with single_qc_file.open("w") as out:
        out.write("sample_id\tfastq\tmedian_len\tmin_len\tmax_len\tqmin\tqmax\tdecision\n")
        for sid, path in sorted(singles.items()):
            q = fastq_qc(path, qc_reads)
            med = q["median_len"] if q["median_len"] != "NA" else 0
            if med >= min_median_len:
                decision = "single_use"
                records.append((sid, path))
            else:
                decision = "skip_too_short"
            out.write(f"{sid}\t{path}\t{q['median_len']}\t{q['min_len']}\t{q['max_len']}\t{q['qmin']}\t{q['qmax']}\t{decision}\n")
    if records:
        analysis_layout = "single_vsearch"
        with manifest_single.open("w") as out:
            out.write("sample-id\tabsolute-filepath\n")
            for sid, path in records:
                out.write(f"{sid}\t{path}\n")
        manifest.write_text(manifest_single.read_text())
        with used_file.open("w") as out:
            out.write("sample_id\tfastq\tmode\n")
            for sid, path in records:
                out.write(f"{sid}\t{path}\tsingle_use\n")
        samples = len(records)
    else:
        analysis_layout = "skip_too_short"
        notes.append("All single-end reads were shorter than MIN_MEDIAN_LEN.")

elif incomplete_pairs and not singles and not complete_pairs:
    input_layout = "incomplete_pairs"
    records = []
    with single_qc_file.open("w") as out:
        out.write("sample_id\tfastq\toriginal_mate\tmedian_len\tmin_len\tmax_len\tqmin\tqmax\tdecision\n")
        for sid, r1, r2 in incomplete_pairs:
            for mate, path in [("R1", r1), ("R2", r2)]:
                if not path:
                    continue
                q = fastq_qc(path, qc_reads)
                med = q["median_len"] if q["median_len"] != "NA" else 0
                if med >= min_median_len:
                    decision = f"single_use_{mate}"
                    records.append((sid, path, decision))
                else:
                    decision = "skip_too_short"
                out.write(f"{sid}\t{path}\t{mate}\t{q['median_len']}\t{q['min_len']}\t{q['max_len']}\t{q['qmin']}\t{q['qmax']}\t{decision}\n")
    if records:
        analysis_layout = "single_best_mate_vsearch"
        with manifest_single.open("w") as out:
            out.write("sample-id\tabsolute-filepath\n")
            for sid, path, decision in records:
                out.write(f"{sid}\t{path}\n")
        manifest.write_text(manifest_single.read_text())
        with used_file.open("w") as out:
            out.write("sample_id\tfastq\tmode\n")
            for sid, path, decision in records:
                out.write(f"{sid}\t{path}\t{decision}\n")
        samples = len(records)
    else:
        analysis_layout = "skip_too_short"
        notes.append("Incomplete pairs found, but no usable long mate.")
else:
    # Mixed project: some samples are single-end and some are paired-end.
    # For cross-PRJNA comparability, convert all usable samples to a unified single-end representation:
    #   - true single files: use the single read if median length passes the threshold;
    #   - paired files with one short index/barcode mate and one long biological mate: use the long mate;
    #   - paired files with two long mates: use the longer mate rather than merging, so the whole mixed project remains one single-end artifact;
    #   - both mates too short: skip that sample.
    input_layout = "mixed"
    analysis_layout = "single_mixed_best_mate_vsearch"
    records = []
    skipped = 0

    with pair_qc_file.open("w") as out:
        out.write("sample_id\tR1\tR2\tR1_median_len\tR2_median_len\tR1_min_len\tR2_min_len\tR1_max_len\tR2_max_len\tR1_qmin\tR1_qmax\tR2_qmin\tR2_qmax\tdecision\n")
        for sid, r1, r2 in complete_pairs:
            q1 = fastq_qc(r1, qc_reads)
            q2 = fastq_qc(r2, qc_reads)
            r1_med = q1["median_len"] if q1["median_len"] != "NA" else 0
            r2_med = q2["median_len"] if q2["median_len"] != "NA" else 0
            r1_long = r1_med >= min_median_len
            r2_long = r2_med >= min_median_len

            if r1_long and r2_long:
                if r2_med > r1_med:
                    decision = "mixed_single_use_R2_both_long_longer"
                    records.append((sid, r2, decision))
                else:
                    decision = "mixed_single_use_R1_both_long_longer"
                    records.append((sid, r1, decision))
            elif r1_long and not r2_long:
                decision = "mixed_single_use_R1"
                records.append((sid, r1, decision))
            elif (not r1_long) and r2_long:
                decision = "mixed_single_use_R2"
                records.append((sid, r2, decision))
            else:
                decision = "mixed_skip_both_short"
                skipped += 1

            out.write(
                f"{sid}\t{r1}\t{r2}\t{q1['median_len']}\t{q2['median_len']}\t"
                f"{q1['min_len']}\t{q2['min_len']}\t{q1['max_len']}\t{q2['max_len']}\t"
                f"{q1['qmin']}\t{q1['qmax']}\t{q2['qmin']}\t{q2['qmax']}\t{decision}\n"
            )

    with single_qc_file.open("w") as out:
        out.write("sample_id\tfastq\tmedian_len\tmin_len\tmax_len\tqmin\tqmax\tdecision\n")
        for sid, path in sorted(singles.items()):
            q = fastq_qc(path, qc_reads)
            med = q["median_len"] if q["median_len"] != "NA" else 0
            if med >= min_median_len:
                decision = "mixed_single_use_original_single"
                records.append((sid, path, decision))
            else:
                decision = "mixed_skip_single_too_short"
                skipped += 1
            out.write(f"{sid}\t{path}\t{q['median_len']}\t{q['min_len']}\t{q['max_len']}\t{q['qmin']}\t{q['qmax']}\t{decision}\n")

    if not records:
        analysis_layout = "skip_too_short"
        notes.append("Mixed project detected, but no usable single or long mate passed MIN_MEDIAN_LEN.")
    else:
        # De-duplicate sample IDs if a project has overlapping single and paired IDs.
        seen = {}
        dedup_records = []
        for sid, path, decision in records:
            if sid in seen:
                seen[sid] += 1
                sid2 = f"{sid}_{seen[sid]}"
            else:
                seen[sid] = 1
                sid2 = sid
            dedup_records.append((sid2, path, decision, sid))

        with manifest_single.open("w") as out:
            out.write("sample-id\tabsolute-filepath\n")
            for sid2, path, decision, original_sid in dedup_records:
                out.write(f"{sid2}\t{path}\n")
        manifest.write_text(manifest_single.read_text())

        with used_file.open("w") as out:
            out.write("sample_id\toriginal_sample_id\tfastq\tmode\n")
            for sid2, path, decision, original_sid in dedup_records:
                out.write(f"{sid2}\t{original_sid}\t{path}\t{decision}\n")

        samples = len(dedup_records)
        notes.append(
            f"Mixed single/paired project converted to unified single-end by selecting each sample's usable biological read. "
            f"usable={samples}, skipped={skipped}."
        )

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
    out.write(f"min_median_len\t{min_median_len}\n")
    out.write(f"search_dir\t{fastq_dir}\n")
    if notes:
        out.write("notes\t" + " ".join(notes).replace("\t", " ") + "\n")

print(layout_file.read_text())
PY

cat "$MANIFEST_DIR/read_layout.tsv" | tee -a "$LOG_DIR/run_qiime2_mate_aware_vsearch.log"

ANALYSIS_LAYOUT="$(awk -F'\t' '$1=="analysis_layout"{print $2}' "$MANIFEST_DIR/read_layout.tsv")"

if [[ "$ANALYSIS_LAYOUT" == "skip_too_short" ]]; then
  echo "[SKIP] $PROJECT has no usable reads with median length >= $MIN_MEDIAN_LEN bp." | tee "$LOG_DIR/skipped.too_short_reads.log"
  exit 0
fi

if [[ "$ANALYSIS_LAYOUT" == "manual_required" ]]; then
  echo "[ERROR] $PROJECT requires manual split/processing. See $MANIFEST_DIR/read_layout.tsv" >&2
  exit 1
fi

if [[ "${AC_STOP_AFTER_QC:-0}" == "1" ]]; then
  log "Read-layout QC completed; classification was not run."
  exit 0
fi

if [[ "$ANALYSIS_LAYOUT" == "paired_merge_then_vsearch" ]]; then
  RAW_DEMUX_QZA="$IMPORT_DIR/paired-end-demux.qza"
  MERGED_DEMUX_QZA="$IMPORT_DIR/merged-demux.qza"
  UNMERGED_DEMUX_QZA="$IMPORT_DIR/unmerged-demux.qza"

  if [[ ! -f "$RAW_DEMUX_QZA" ]]; then
    log "QIIME2 import as PairedEndFastqManifestPhred33V2"
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path "$MANIFEST_DIR/manifest.tsv" \
      --output-path "$RAW_DEMUX_QZA" \
      --input-format PairedEndFastqManifestPhred33V2
  else
    log "Paired import output exists; skip"
  fi

  if [[ ! -f "$DEMUX_DIR/demux-raw-summary.qzv" ]]; then
    log "demux summarize: raw paired-end"
    qiime demux summarize \
      --i-data "$RAW_DEMUX_QZA" \
      --o-visualization "$DEMUX_DIR/demux-raw-summary.qzv"
  fi

  if [[ ! -f "$MERGED_DEMUX_QZA" ]]; then
    log "VSEARCH merge-pairs for true paired-end reads"
    qiime vsearch merge-pairs \
      --i-demultiplexed-seqs "$RAW_DEMUX_QZA" \
      --p-threads "$THREADS" \
      --o-merged-sequences "$MERGED_DEMUX_QZA" \
      --o-unmerged-sequences "$UNMERGED_DEMUX_QZA"
  else
    log "Merged demux exists; skip"
  fi

  if [[ ! -f "$DEMUX_DIR/demux-merged-summary.qzv" ]]; then
    log "demux summarize: merged paired-end"
    qiime demux summarize \
      --i-data "$MERGED_DEMUX_QZA" \
      --o-visualization "$DEMUX_DIR/demux-merged-summary.qzv"
  fi

  DEMUX_FOR_VSEARCH="$MERGED_DEMUX_QZA"

else
  RAW_DEMUX_QZA="$IMPORT_DIR/single-end-demux.qza"

  if [[ ! -f "$RAW_DEMUX_QZA" ]]; then
    log "QIIME2 import as SingleEndFastqManifestPhred33V2"
    qiime tools import \
      --type 'SampleData[SequencesWithQuality]' \
      --input-path "$MANIFEST_DIR/manifest.tsv" \
      --output-path "$RAW_DEMUX_QZA" \
      --input-format SingleEndFastqManifestPhred33V2
  else
    log "Single import output exists; skip"
  fi

  if [[ ! -f "$DEMUX_DIR/demux-raw-summary.qzv" ]]; then
    log "demux summarize: single-end / selected long mate"
    qiime demux summarize \
      --i-data "$RAW_DEMUX_QZA" \
      --o-visualization "$DEMUX_DIR/demux-raw-summary.qzv"
  fi

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

if [[ ! -f "$DEREP_TABLE" || ! -f "$DEREP_SEQS" ]]; then
  log "VSEARCH dereplicate sequences"
  qiime vsearch dereplicate-sequences \
    --i-sequences "$DEMUX_FOR_VSEARCH" \
    --o-dereplicated-table "$DEREP_TABLE" \
    --o-dereplicated-sequences "$DEREP_SEQS"
else
  log "VSEARCH dereplicated outputs exist; skip"
fi

if [[ ! -f "$CLUSTER_TABLE" || ! -f "$CLUSTER_SEQS" ]]; then
  log "VSEARCH de novo clustering at identity=$PERC_IDENTITY"
  qiime vsearch cluster-features-de-novo \
    --i-table "$DEREP_TABLE" \
    --i-sequences "$DEREP_SEQS" \
    --p-perc-identity "$PERC_IDENTITY" \
    --p-threads "$THREADS" \
    --o-clustered-table "$CLUSTER_TABLE" \
    --o-clustered-sequences "$CLUSTER_SEQS"
else
  log "VSEARCH clustered outputs exist; skip"
fi

if [[ ! -f "$NONCHIMERAS" ]]; then
  log "VSEARCH uchime-denovo chimera detection"
  qiime vsearch uchime-denovo \
    --i-table "$CLUSTER_TABLE" \
    --i-sequences "$CLUSTER_SEQS" \
    --o-chimeras "$CHIMERAS" \
    --o-nonchimeras "$NONCHIMERAS" \
    --o-stats "$UCHIME_STATS"
else
  log "UCHIME outputs exist; skip"
fi

if [[ ! -f "$TABLE_QZA" ]]; then
  log "Filter chimeric features from table"
  qiime feature-table filter-features \
    --i-table "$CLUSTER_TABLE" \
    --m-metadata-file "$NONCHIMERAS" \
    --o-filtered-table "$TABLE_QZA"
else
  log "Final table exists; skip"
fi

if [[ ! -f "$REP_SEQS_QZA" ]]; then
  log "Filter chimeric features from representative sequences"
  qiime feature-table filter-seqs \
    --i-data "$CLUSTER_SEQS" \
    --m-metadata-file "$NONCHIMERAS" \
    --o-filtered-data "$REP_SEQS_QZA"
else
  log "Final rep-seqs exists; skip"
fi

log "Generate feature table and rep-seqs summaries"
qiime feature-table summarize \
  --i-table "$TABLE_QZA" \
  --o-visualization "$VSEARCH_DIR/table.qzv" || true

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

cat <<EOF2
============================================================
Finished mate-aware QIIME2 VSEARCH workflow: $PROJECT

Main outputs:
  $EXPORT_DIR/feature-table.tsv
  $TAX_DIR/taxonomy.tsv
  $VSEARCH_DIR/table.qza
  $VSEARCH_DIR/rep-seqs.qza

Layout/QC:
  $MANIFEST_DIR/read_layout.tsv
  $MANIFEST_DIR/fastq_used.tsv
  $MANIFEST_DIR/pair_mate_length_qc.tsv
  $MANIFEST_DIR/single_length_qc.tsv
============================================================
EOF2
