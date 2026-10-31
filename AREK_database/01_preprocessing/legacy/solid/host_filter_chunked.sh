#!/usr/bin/env bash
set -Eeuo pipefail
# Portable recovered source. Canonical resumable route: run_pipeline.py.
: "${RAW_DATA_ROOT:?set raw data root}" "${RESULT_ROOT:?set isolated output root}" "${BFAST_HOME:?set installed BFAST 0.5.6 prefix}" "${SOLID_ENV:?set dedicated colorspace environment prefix}" "${BFAST_DB_ROOT:?set BFAST hg38 CS database root}"
if [[ "${AREK_LEGACY_SOURCE_MODE:-}" != "isolated_research" ]]; then echo "Recovered historical entry; use canonical run_pipeline.py. To inspect native batching/diagnostics in a NEW output directory set AREK_LEGACY_SOURCE_MODE=isolated_research." >&2; exit 2; fi


trap '
rc=$?
echo "[ERROR] Failed at line ${LINENO}: ${BASH_COMMAND}" >&2
echo "[ERROR] Exit code: ${rc}" >&2
exit "${rc}"
' ERR

###############################################################################
# Basic configuration
###############################################################################

unset LD_LIBRARY_PATH

BFAST_HOME="${BFAST_HOME}"
SOLID_ENV="${SOLID_ENV}"

BMATCHES="${BFAST_HOME}/bin/bmatches"
BALIGN="${BFAST_HOME}/bin/balign"
BPOSTPROCESS="${BFAST_HOME}/bin/bpostprocess"

PYTHON="${SOLID_ENV}/bin/python"
SAMTOOLS="${SOLID_ENV}/bin/samtools"
TIME_BIN="/usr/bin/time"

PROJECT="PRJEB18041"

SAMPLE="${1:-${SAMPLE:-ADS_10}}"

RAW_FASTQ="${RAW_FASTQ:-${RAW_DATA_ROOT}/${PROJECT}/fastq/${SAMPLE}_r.fastq.gz}"

# 0 means use the complete sample.
N_READS="${N_READS:-0}"

# Number of reads processed by one BFAST invocation.
CHUNK_READS="${CHUNK_READS:-100000}"

# BFAST 0.5.6 is most stable with one thread on this server.
THREADS="${THREADS:-1}"

# FORCE=1 deletes the corresponding run and starts again.
FORCE="${FORCE:-0}"

# Delete BMF, BAF, SAM and prepared chunk after that chunk passes QC.
CLEAN_INTERMEDIATE="${CLEAN_INTERMEDIATE:-1}"

OUTPUT_BASE="${OUTPUT_BASE:-${RESULT_ROOT}/${PROJECT}/host_filter_chunked}"

###############################################################################
# Validate numeric parameters
###############################################################################

for value_name in N_READS CHUNK_READS THREADS FORCE CLEAN_INTERMEDIATE; do
    value="${!value_name}"

    if [[ ! "${value}" =~ ^[0-9]+$ ]]; then
        echo "[ERROR] ${value_name} must be a non-negative integer: ${value}" >&2
        exit 1
    fi
done

if [[ "${CHUNK_READS}" -le 0 ]]; then
    echo "[ERROR] CHUNK_READS must be greater than zero" >&2
    exit 1
fi

if [[ "${THREADS}" -le 0 ]]; then
    echo "[ERROR] THREADS must be greater than zero" >&2
    exit 1
fi

###############################################################################
# hg38 BFAST index
###############################################################################

INDEX_DIR="${BFAST_DB_ROOT}/index"

NT_BRG="${INDEX_DIR}/bfast.rg.file.hg38.0.brg"
CS_BRG="${INDEX_DIR}/bfast.rg.file.hg38.1.brg"
MAIN_INDEXES="${INDEX_DIR}/main.indexes.hg38.txt"
SECONDARY_INDEXES="${INDEX_DIR}/secondary.indexes.hg38.txt"
OFFSETS="${INDEX_DIR}/offsets.0-100.txt"
SCORING="${INDEX_DIR}/scoring.1.txt"

if [[ "${N_READS}" -eq 0 ]]; then
    READ_TAG="all"
else
    READ_TAG="n${N_READS}"
fi

RUN_ID="${SAMPLE}.hg38.${READ_TAG}.chunk${CHUNK_READS}.qfix"

RUN_ROOT="${OUTPUT_BASE}/${RUN_ID}"
CHUNK_DIR="${RUN_ROOT}/prepared_chunks"
CHUNK_RESULT_DIR="${RUN_ROOT}/chunk_results"
LOG_DIR="${RUN_ROOT}/logs"
FINAL_DIR="${RUN_ROOT}/final"

MANIFEST="${RUN_ROOT}/chunk_manifest.tsv"
PREPARE_SUMMARY="${RUN_ROOT}/prepare_summary.tsv"
PREPARE_DONE="${RUN_ROOT}/.prepare.done"

FINAL_FASTQ="${FINAL_DIR}/${RUN_ID}.nonhuman.bfast.colorspace.fastq.gz"
FINAL_MAPPED_IDS="${FINAL_DIR}/${RUN_ID}.human_mapped_read_ids.txt"
FINAL_SUMMARY="${FINAL_DIR}/${RUN_ID}.host_filter_summary.tsv"

###############################################################################
# Initialization
###############################################################################

echo "[INFO] Start: $(date)"
echo "[INFO] PROJECT=${PROJECT}"
echo "[INFO] SAMPLE=${SAMPLE}"
echo "[INFO] RAW_FASTQ=${RAW_FASTQ}"
echo "[INFO] N_READS=${N_READS}"
echo "[INFO] CHUNK_READS=${CHUNK_READS}"
echo "[INFO] THREADS=${THREADS}"
echo "[INFO] FORCE=${FORCE}"
echo "[INFO] CLEAN_INTERMEDIATE=${CLEAN_INTERMEDIATE}"
echo "[INFO] RUN_ID=${RUN_ID}"
echo "[INFO] RUN_ROOT=${RUN_ROOT}"

if [[ "${FORCE}" == "1" ]]; then
    echo "[INFO] FORCE=1: removing previous run"
    rm -rf "${RUN_ROOT}"
fi

mkdir -p \
    "${OUTPUT_BASE}" \
    "${RUN_ROOT}" \
    "${CHUNK_DIR}" \
    "${CHUNK_RESULT_DIR}" \
    "${LOG_DIR}" \
    "${FINAL_DIR}"

if [[ -s "${FINAL_SUMMARY}" && "${FORCE}" == "0" ]]; then
    if grep -q $'accounting_status\tPASS' "${FINAL_SUMMARY}"; then
        echo "[SKIP] Complete PASS result already exists:"
        echo "${FINAL_SUMMARY}"
        exit 0
    fi
fi

for executable in \
    "${BMATCHES}" \
    "${BALIGN}" \
    "${BPOSTPROCESS}" \
    "${PYTHON}" \
    "${SAMTOOLS}" \
    "${TIME_BIN}"
do
    if [[ ! -x "${executable}" ]]; then
        echo "[ERROR] Executable not found: ${executable}" >&2
        exit 1
    fi
done

for input_file in \
    "${RAW_FASTQ}" \
    "${NT_BRG}" \
    "${CS_BRG}" \
    "${MAIN_INDEXES}" \
    "${SECONDARY_INDEXES}" \
    "${OFFSETS}" \
    "${SCORING}"
do
    if [[ ! -s "${input_file}" ]]; then
        echo "[ERROR] Missing or empty input: ${input_file}" >&2
        exit 1
    fi
done

###############################################################################
# Step 0: Split and prepare color-space FASTQ
#
# Original:
#   sequence = primer base + 50 colors, length 51
#   quality  = length 51
#
# BFAST-compatible:
#   sequence length 51
#   quality length 50
#
# "." is changed to "4".
# The first quality character is removed when QUAL length == SEQ length.
###############################################################################

if [[ -s "${MANIFEST}" && -f "${PREPARE_DONE}" && "${FORCE}" == "0" ]]; then
    echo "[SKIP] Prepared chunk manifest already exists"
else
    echo "[INFO] Step 0: preparing color-space chunks"

    rm -rf "${CHUNK_DIR}" "${CHUNK_RESULT_DIR}"
    mkdir -p "${CHUNK_DIR}" "${CHUNK_RESULT_DIR}"

    rm -f \
        "${MANIFEST}" \
        "${PREPARE_SUMMARY}" \
        "${PREPARE_DONE}"

    "${PYTHON}" - \
        "${RAW_FASTQ}" \
        "${CHUNK_DIR}" \
        "${MANIFEST}" \
        "${PREPARE_SUMMARY}" \
        "${CHUNK_READS}" \
        "${N_READS}" <<'PY'
import gzip
import re
import sys
from collections import Counter
from pathlib import Path

source = Path(sys.argv[1])
chunk_dir = Path(sys.argv[2])
manifest_path = Path(sys.argv[3])
summary_path = Path(sys.argv[4])
chunk_size = int(sys.argv[5])
read_limit = int(sys.argv[6])

valid_sequence = re.compile(r"^[ACGTNacgtn][0-4]+$")

chunk_dir.mkdir(parents=True, exist_ok=True)

seen = 0
written = 0
bad_sequence = 0
bad_length = 0
quality_first_character_removed = 0

length_patterns = Counter()
manifest_rows = []

chunk_index = 0
chunk_count = 0
chunk_handle = None
chunk_path = None


def close_current_chunk():
    global chunk_handle, chunk_path, chunk_count

    if chunk_handle is not None:
        chunk_handle.close()

        manifest_rows.append(
            (
                chunk_index,
                f"chunk_{chunk_index:06d}",
                chunk_count,
                str(chunk_path),
            )
        )

        chunk_handle = None
        chunk_path = None
        chunk_count = 0


with gzip.open(source, "rt", errors="replace") as fin:
    while True:
        if read_limit > 0 and seen >= read_limit:
            break

        header = fin.readline()

        if not header:
            break

        sequence = fin.readline().rstrip("\r\n")
        plus = fin.readline()
        raw_quality = fin.readline().rstrip("\r\n")

        if not plus or raw_quality == "":
            raise RuntimeError(
                f"Truncated FASTQ record after read {seen}"
            )

        seen += 1
        sequence = sequence.replace(".", "4")

        if not valid_sequence.fullmatch(sequence):
            bad_sequence += 1
            continue

        if len(raw_quality) == len(sequence):
            quality = raw_quality[1:]
            quality_first_character_removed += 1
        elif len(raw_quality) == len(sequence) - 1:
            quality = raw_quality
        else:
            bad_length += 1
            continue

        if len(quality) != len(sequence) - 1:
            bad_length += 1
            continue

        if chunk_handle is None:
            chunk_index += 1
            chunk_path = (
                chunk_dir /
                f"chunk_{chunk_index:06d}.bfast.colorspace.fastq"
            )
            chunk_handle = chunk_path.open("w")
            chunk_count = 0

        chunk_handle.write(header)
        chunk_handle.write(sequence + "\n")
        chunk_handle.write(plus)
        chunk_handle.write(quality + "\n")

        chunk_count += 1
        written += 1

        length_patterns[
            (len(sequence), len(raw_quality), len(quality))
        ] += 1

        if chunk_count >= chunk_size:
            close_current_chunk()

close_current_chunk()

if written == 0:
    raise SystemExit("No valid BFAST-compatible reads were generated")

with manifest_path.open("w") as manifest:
    manifest.write("chunk_number\tchunk_id\treads\tfastq\n")

    for row in manifest_rows:
        manifest.write(
            f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\n"
        )

with summary_path.open("w") as summary:
    summary.write("metric\tvalue\n")
    summary.write(f"reads_seen\t{seen}\n")
    summary.write(f"reads_written\t{written}\n")
    summary.write(f"bad_sequence\t{bad_sequence}\n")
    summary.write(f"bad_length\t{bad_length}\n")
    summary.write(
        "quality_first_character_removed\t"
        f"{quality_first_character_removed}\n"
    )
    summary.write(f"chunks\t{len(manifest_rows)}\n")
    summary.write(f"chunk_reads\t{chunk_size}\n")

    for pattern, count in sorted(length_patterns.items()):
        seq_len, raw_q_len, bfast_q_len = pattern
        summary.write(
            "length_pattern\t"
            f"sequence:{seq_len},"
            f"raw_quality:{raw_q_len},"
            f"bfast_quality:{bfast_q_len},"
            f"reads:{count}\n"
        )

print(f"reads_seen={seen}")
print(f"reads_written={written}")
print(f"bad_sequence={bad_sequence}")
print(f"bad_length={bad_length}")
print(f"chunks={len(manifest_rows)}")
print(f"chunk_reads={chunk_size}")
PY

    touch "${PREPARE_DONE}"
fi

TOTAL_CHUNKS=$(awk 'NR > 1 {n++} END {print n+0}' "${MANIFEST}")
TOTAL_PREPARED=$(awk -F '\t' 'NR > 1 {n += $3} END {print n+0}' "${MANIFEST}")

echo "[INFO] Prepared reads: ${TOTAL_PREPARED}"
echo "[INFO] Total chunks: ${TOTAL_CHUNKS}"

###############################################################################
# Process one chunk at a time
###############################################################################

while IFS=$'\t' read -r chunk_number chunk_id chunk_reads chunk_fastq
do
    if [[ "${chunk_number}" == "chunk_number" ]]; then
        continue
    fi

    BFAST_ID="${SAMPLE}.${chunk_id}.hg38.qfix"

    CHUNK_ROOT="${CHUNK_RESULT_DIR}/${chunk_id}"
    OUT_DIR="${CHUNK_ROOT}/output"
    SAM_DIR="${CHUNK_ROOT}/sam"
    TMP_DIR="${CHUNK_ROOT}/tmp"
    CHUNK_LOG_DIR="${CHUNK_ROOT}/logs"

    MATCH_FILE="${OUT_DIR}/bfast.matches.file.${BFAST_ID}.bmf"
    ALIGN_FILE="${OUT_DIR}/bfast.aligned.file.${BFAST_ID}.baf"
    SAM_FILE="${SAM_DIR}/bfast.reported.file.${BFAST_ID}.sam"

    MAPPED_IDS="${OUT_DIR}/${BFAST_ID}.human_mapped_read_ids.txt"
    NONHUMAN_FASTQ="${OUT_DIR}/${BFAST_ID}.nonhuman.colorspace.fastq.gz"
    CHUNK_SUMMARY="${OUT_DIR}/${BFAST_ID}.host_filter_summary.tsv"
    UNRESOLVED_IDS="${OUT_DIR}/${BFAST_ID}.unresolved_mapped_ids.txt"
    DONE_FILE="${CHUNK_ROOT}/.done"

    if [[ -f "${DONE_FILE}" && -s "${CHUNK_SUMMARY}" ]]; then
        if grep -q $'accounting_status\tPASS' "${CHUNK_SUMMARY}"; then
            echo "[SKIP] ${chunk_id}: already completed"
            continue
        fi
    fi

    if [[ ! -s "${chunk_fastq}" ]]; then
        echo "[ERROR] Missing prepared chunk: ${chunk_fastq}" >&2
        exit 1
    fi

    mkdir -p \
        "${OUT_DIR}" \
        "${SAM_DIR}" \
        "${TMP_DIR}" \
        "${CHUNK_LOG_DIR}"

    echo
    echo "[INFO] Processing ${chunk_id}/${TOTAL_CHUNKS}"
    echo "[INFO] Reads=${chunk_reads}"
    echo "[INFO] Started=$(date)"

    ###########################################################################
    # bmatches
    ###########################################################################

    rm -f "${MATCH_FILE}"

    "${TIME_BIN}" -v \
    "${BMATCHES}" \
        -r "${CS_BRG}" \
        -i "${MAIN_INDEXES}" \
        -I "${SECONDARY_INDEXES}" \
        -O "${OFFSETS}" \
        -R "${chunk_fastq}" \
        -A 1 \
        -n "${THREADS}" \
        -o "${BFAST_ID}" \
        -d "${OUT_DIR}/" \
        -T "${TMP_DIR}/" \
        > "${CHUNK_LOG_DIR}/01_bmatches.stdout.log" \
        2> "${CHUNK_LOG_DIR}/01_bmatches.stderr.log"

    if [[ ! -s "${MATCH_FILE}" ]]; then
        echo "[ERROR] Missing bmatches output for ${chunk_id}" >&2
        exit 1
    fi

    ###########################################################################
    # balign
    ###########################################################################

    rm -f \
        "${ALIGN_FILE}" \
        "${OUT_DIR}/bfast.not.aligned.file.${BFAST_ID}.bmf"

    "${TIME_BIN}" -v \
    "${BALIGN}" \
        -r "${NT_BRG}" \
        -m "${MATCH_FILE}" \
        -x "${SCORING}" \
        -A 1 \
        -O 15 \
        -n "${THREADS}" \
        -o "${BFAST_ID}" \
        -d "${OUT_DIR}/" \
        -T "${TMP_DIR}/" \
        > "${CHUNK_LOG_DIR}/02_balign.stdout.log" \
        2> "${CHUNK_LOG_DIR}/02_balign.stderr.log"

    if [[ ! -s "${ALIGN_FILE}" ]]; then
        echo "[ERROR] Missing balign output for ${chunk_id}" >&2
        exit 1
    fi

    ###########################################################################
    # bpostprocess
    #
    # -a 1 retained:
    # -a 3 caused a segmentation fault in this BFAST build.
    #
    # -O 3 writes SAM directly.
    ###########################################################################

    rm -f \
        "${SAM_FILE}" \
        "${SAM_DIR}/bfast.not.reported.file.${BFAST_ID}.baf"

    "${TIME_BIN}" -v \
    "${BPOSTPROCESS}" \
        -r "${NT_BRG}" \
        -i "${ALIGN_FILE}" \
        -a 1 \
        -O 3 \
        -o "${BFAST_ID}" \
        -d "${SAM_DIR}/" \
        > "${CHUNK_LOG_DIR}/03_bpostprocess.stdout.log" \
        2> "${CHUNK_LOG_DIR}/03_bpostprocess.stderr.log"

    if [[ ! -s "${SAM_FILE}" ]]; then
        echo "[ERROR] Missing SAM output for ${chunk_id}" >&2
        exit 1
    fi

    if ! "${SAMTOOLS}" view -c "${SAM_FILE}" \
        > "${CHUNK_LOG_DIR}/03_sam_record_count.txt" \
        2> "${CHUNK_LOG_DIR}/03_sam_validation.stderr.log"
    then
        echo "[ERROR] Invalid SAM for ${chunk_id}" >&2
        echo "[ERROR] See ${CHUNK_LOG_DIR}/03_sam_validation.stderr.log" >&2
        exit 1
    fi

    ###########################################################################
    # Normalize BFAST QNAMEs and filter host reads
    ###########################################################################

    rm -f \
        "${MAPPED_IDS}" \
        "${NONHUMAN_FASTQ}" \
        "${CHUNK_SUMMARY}" \
        "${UNRESOLVED_IDS}"

    "${PYTHON}" - \
        "${BFAST_ID}" \
        "${chunk_id}" \
        "${chunk_fastq}" \
        "${SAM_FILE}" \
        "${MAPPED_IDS}" \
        "${NONHUMAN_FASTQ}" \
        "${CHUNK_SUMMARY}" \
        "${UNRESOLVED_IDS}" <<'PY'
import gzip
import sys
from pathlib import Path

prefix = sys.argv[1]
chunk_id = sys.argv[2]

fastq_path = Path(sys.argv[3])
sam_path = Path(sys.argv[4])
mapped_path = Path(sys.argv[5])
nonhuman_path = Path(sys.argv[6])
summary_path = Path(sys.argv[7])
unresolved_path = Path(sys.argv[8])


def remove_mate_suffix(name):
    if name.endswith("/1") or name.endswith("/2"):
        return name[:-2]
    return name


def fastq_records(path):
    with path.open() as handle:
        while True:
            header = handle.readline()

            if not header:
                break

            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not quality:
                raise RuntimeError(
                    f"Truncated FASTQ record in {path}"
                )

            yield header, sequence, plus, quality


fastq_ids = set()

for header, _, _, _ in fastq_records(fastq_path):
    read_id = header[1:].strip().split()[0]
    fastq_ids.add(remove_mate_suffix(read_id))


def normalize_qname(name):
    name = remove_mate_suffix(name.strip())

    candidates = [name]

    if name.startswith(prefix):
        remainder = name[len(prefix):]
        remainder = remainder.lstrip(":._-/|")
        candidates.append(remove_mate_suffix(remainder))

    for separator in (":", "|"):
        if separator in name:
            candidates.append(
                remove_mate_suffix(name.rsplit(separator, 1)[-1])
            )

    for candidate in candidates:
        if candidate in fastq_ids:
            return candidate

    return None


sam_records = 0
mapped_sam_records = 0
mapped_ids = set()
unresolved_names = set()

with sam_path.open(errors="replace") as sam_handle:
    for line_number, line in enumerate(sam_handle, start=1):
        if not line or line.startswith("@"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 11:
            raise RuntimeError(
                f"Malformed SAM line {line_number}"
            )

        try:
            flag = int(fields[1])
        except ValueError as exc:
            raise RuntimeError(
                f"Invalid SAM flag at line {line_number}"
            ) from exc

        sam_records += 1

        if flag & 4:
            continue

        mapped_sam_records += 1

        normalized = normalize_qname(fields[0])

        if normalized is None:
            unresolved_names.add(fields[0])
        else:
            mapped_ids.add(normalized)

if unresolved_names:
    with unresolved_path.open("w") as handle:
        for name in sorted(unresolved_names):
            handle.write(name + "\n")

    raise RuntimeError(
        f"{len(unresolved_names)} mapped SAM QNAMEs could not "
        f"be matched to FASTQ IDs; see {unresolved_path}"
    )

with mapped_path.open("w") as handle:
    for read_id in sorted(mapped_ids):
        handle.write(read_id + "\n")

total_reads = 0
human_removed = 0
nonhuman_retained = 0
found_mapped_ids = set()

with gzip.open(nonhuman_path, "wt") as output_handle:
    for header, sequence, plus, quality in fastq_records(fastq_path):
        total_reads += 1

        read_id = header[1:].strip().split()[0]
        read_id = remove_mate_suffix(read_id)

        if read_id in mapped_ids:
            human_removed += 1
            found_mapped_ids.add(read_id)
            continue

        output_handle.write(header)
        output_handle.write(sequence)
        output_handle.write(plus)
        output_handle.write(quality)

        nonhuman_retained += 1

missing_ids = mapped_ids - found_mapped_ids

if missing_ids:
    with unresolved_path.open("w") as handle:
        for read_id in sorted(missing_ids):
            handle.write(read_id + "\n")

    raise RuntimeError(
        f"{len(missing_ids)} normalized mapped IDs were not found "
        f"during FASTQ filtering"
    )

unique_human_mapped_reads = len(mapped_ids)

if unique_human_mapped_reads != human_removed:
    raise RuntimeError(
        "Mapped-read accounting failed: "
        f"mapped={unique_human_mapped_reads}, "
        f"removed={human_removed}"
    )

if total_reads != human_removed + nonhuman_retained:
    raise RuntimeError(
        "FASTQ accounting failed: "
        f"input={total_reads}, "
        f"removed={human_removed}, "
        f"retained={nonhuman_retained}"
    )

human_fraction = (
    human_removed / total_reads
    if total_reads
    else 0.0
)

with summary_path.open("w") as summary:
    summary.write("metric\tvalue\n")
    summary.write(f"chunk_id\t{chunk_id}\n")
    summary.write(f"input_reads\t{total_reads}\n")
    summary.write(f"sam_alignment_records\t{sam_records}\n")
    summary.write(f"mapped_sam_records\t{mapped_sam_records}\n")
    summary.write(
        f"unique_human_mapped_reads\t"
        f"{unique_human_mapped_reads}\n"
    )
    summary.write(f"human_removed\t{human_removed}\n")
    summary.write(
        f"human_mapped_fraction\t{human_fraction:.8f}\n"
    )
    summary.write(
        f"nonhuman_retained_reads\t{nonhuman_retained}\n"
    )
    summary.write("accounting_status\tPASS\n")

print(f"chunk_id={chunk_id}")
print(f"input_reads={total_reads}")
print(
    f"unique_human_mapped_reads="
    f"{unique_human_mapped_reads}"
)
print(f"human_removed={human_removed}")
print(f"nonhuman_retained={nonhuman_retained}")
print("accounting_status=PASS")
PY

    if ! grep -q $'accounting_status\tPASS' "${CHUNK_SUMMARY}"; then
        echo "[ERROR] Accounting did not pass for ${chunk_id}" >&2
        exit 1
    fi

    gzip -t "${NONHUMAN_FASTQ}"

    ###########################################################################
    # Cleanup after successful chunk
    ###########################################################################

    if [[ "${CLEAN_INTERMEDIATE}" == "1" ]]; then
        rm -f \
            "${MATCH_FILE}" \
            "${ALIGN_FILE}" \
            "${SAM_FILE}" \
            "${OUT_DIR}/bfast.not.aligned.file.${BFAST_ID}.bmf" \
            "${SAM_DIR}/bfast.not.reported.file.${BFAST_ID}.baf" \
            "${chunk_fastq}"

        rm -rf "${TMP_DIR}"
    fi

    touch "${DONE_FILE}"

    echo "[INFO] Completed ${chunk_id}: $(date)"

done < "${MANIFEST}"

###############################################################################
# Verify all chunks completed
###############################################################################

DONE_CHUNKS=0

while IFS=$'\t' read -r chunk_number chunk_id chunk_reads chunk_fastq
do
    if [[ "${chunk_number}" == "chunk_number" ]]; then
        continue
    fi

    if [[ -f "${CHUNK_RESULT_DIR}/${chunk_id}/.done" ]]; then
        DONE_CHUNKS=$((DONE_CHUNKS + 1))
    fi
done < "${MANIFEST}"

if [[ "${DONE_CHUNKS}" -ne "${TOTAL_CHUNKS}" ]]; then
    echo "[ERROR] Not all chunks completed:" >&2
    echo "[ERROR] completed=${DONE_CHUNKS}, total=${TOTAL_CHUNKS}" >&2
    exit 1
fi

###############################################################################
# Merge nonhuman FASTQ and mapped IDs
###############################################################################

echo
echo "[INFO] All chunks completed"
echo "[INFO] Merging final outputs"

FINAL_FASTQ_TMP="${FINAL_FASTQ}.tmp.$$"
FINAL_IDS_TMP="${FINAL_MAPPED_IDS}.tmp.$$"

rm -f "${FINAL_FASTQ_TMP}" "${FINAL_IDS_TMP}"
: > "${FINAL_FASTQ_TMP}"
: > "${FINAL_IDS_TMP}"

while IFS=$'\t' read -r chunk_number chunk_id chunk_reads chunk_fastq
do
    if [[ "${chunk_number}" == "chunk_number" ]]; then
        continue
    fi

    BFAST_ID="${SAMPLE}.${chunk_id}.hg38.qfix"
    CHUNK_OUT="${CHUNK_RESULT_DIR}/${chunk_id}/output"

    CHUNK_NONHUMAN="${CHUNK_OUT}/${BFAST_ID}.nonhuman.colorspace.fastq.gz"
    CHUNK_IDS="${CHUNK_OUT}/${BFAST_ID}.human_mapped_read_ids.txt"

    if [[ ! -s "${CHUNK_NONHUMAN}" ]]; then
        echo "[ERROR] Missing chunk FASTQ: ${CHUNK_NONHUMAN}" >&2
        exit 1
    fi

    cat "${CHUNK_NONHUMAN}" >> "${FINAL_FASTQ_TMP}"

    if [[ -s "${CHUNK_IDS}" ]]; then
        cat "${CHUNK_IDS}" >> "${FINAL_IDS_TMP}"
    fi

done < "${MANIFEST}"

gzip -t "${FINAL_FASTQ_TMP}"

mv "${FINAL_FASTQ_TMP}" "${FINAL_FASTQ}"

LC_ALL=C sort -u "${FINAL_IDS_TMP}" > "${FINAL_MAPPED_IDS}"
rm -f "${FINAL_IDS_TMP}"

###############################################################################
# Aggregate final summary
###############################################################################

"${PYTHON}" - \
    "${MANIFEST}" \
    "${CHUNK_RESULT_DIR}" \
    "${SAMPLE}" \
    "${N_READS}" \
    "${CHUNK_READS}" \
    "${THREADS}" \
    "${FINAL_MAPPED_IDS}" \
    "${FINAL_SUMMARY}" <<'PY'
import sys
from pathlib import Path

manifest_path = Path(sys.argv[1])
chunk_result_dir = Path(sys.argv[2])
sample = sys.argv[3]
requested_reads = int(sys.argv[4])
chunk_reads = int(sys.argv[5])
threads = int(sys.argv[6])
final_mapped_ids = Path(sys.argv[7])
final_summary = Path(sys.argv[8])

total_input = 0
total_sam_records = 0
total_mapped_records = 0
total_removed = 0
total_retained = 0
total_chunks = 0

with manifest_path.open() as manifest:
    header = manifest.readline()

    for line in manifest:
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 4:
            continue

        chunk_id = fields[1]
        bfast_id = f"{sample}.{chunk_id}.hg38.qfix"

        summary_path = (
            chunk_result_dir /
            chunk_id /
            "output" /
            f"{bfast_id}.host_filter_summary.tsv"
        )

        if not summary_path.is_file():
            raise RuntimeError(
                f"Missing chunk summary: {summary_path}"
            )

        metrics = {}

        with summary_path.open() as summary:
            next(summary)

            for row in summary:
                key, value = row.rstrip("\n").split("\t", 1)
                metrics[key] = value

        if metrics.get("accounting_status") != "PASS":
            raise RuntimeError(
                f"Chunk did not pass accounting: {chunk_id}"
            )

        total_chunks += 1
        total_input += int(metrics["input_reads"])
        total_sam_records += int(metrics["sam_alignment_records"])
        total_mapped_records += int(metrics["mapped_sam_records"])
        total_removed += int(metrics["human_removed"])
        total_retained += int(metrics["nonhuman_retained_reads"])

with final_mapped_ids.open() as handle:
    unique_final_mapped = sum(
        1 for line in handle if line.strip()
    )

if unique_final_mapped != total_removed:
    raise RuntimeError(
        "Final mapped-ID accounting failed: "
        f"unique_final_mapped={unique_final_mapped}, "
        f"removed={total_removed}"
    )

if total_input != total_removed + total_retained:
    raise RuntimeError(
        "Final FASTQ accounting failed: "
        f"input={total_input}, "
        f"removed={total_removed}, "
        f"retained={total_retained}"
    )

fraction = total_removed / total_input if total_input else 0.0

with final_summary.open("w") as summary:
    summary.write("metric\tvalue\n")
    summary.write(f"sample\t{sample}\n")
    summary.write(f"requested_reads\t{requested_reads}\n")
    summary.write(f"chunk_reads\t{chunk_reads}\n")
    summary.write(f"completed_chunks\t{total_chunks}\n")
    summary.write(f"input_reads\t{total_input}\n")
    summary.write(
        f"sam_alignment_records\t{total_sam_records}\n"
    )
    summary.write(
        f"mapped_sam_records\t{total_mapped_records}\n"
    )
    summary.write(
        f"unique_human_mapped_reads\t"
        f"{unique_final_mapped}\n"
    )
    summary.write(f"human_removed\t{total_removed}\n")
    summary.write(
        f"human_mapped_fraction\t{fraction:.8f}\n"
    )
    summary.write(
        f"nonhuman_retained_reads\t{total_retained}\n"
    )
    summary.write(f"threads\t{threads}\n")
    summary.write("accounting_status\tPASS\n")

print(f"input_reads={total_input}")
print(f"human_removed={total_removed}")
print(f"nonhuman_retained={total_retained}")
print(f"human_mapped_fraction={fraction:.8f}")
print("accounting_status=PASS")
PY

echo
echo "[INFO] Completed successfully: $(date)"
echo

if command -v column >/dev/null 2>&1; then
    column -t -s $'\t' "${FINAL_SUMMARY}"
else
    cat "${FINAL_SUMMARY}"
fi

echo
echo "[INFO] Final files:"
ls -lh \
    "${FINAL_FASTQ}" \
    "${FINAL_MAPPED_IDS}" \
    "${FINAL_SUMMARY}"

echo "[INFO] Done."
