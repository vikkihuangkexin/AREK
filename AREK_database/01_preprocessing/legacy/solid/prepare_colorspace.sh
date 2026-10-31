#!/usr/bin/env bash
set -euo pipefail
# Portable recovered source. Canonical resumable route: run_pipeline.py.
: "${RAW_DATA_ROOT:?set raw data root}" "${RESULT_ROOT:?set isolated output root}" "${BFAST_HOME:?set installed BFAST 0.5.6 prefix}" "${SOLID_ENV:?set dedicated colorspace environment prefix}" "${BFAST_DB_ROOT:?set BFAST hg38 CS database root}"
if [[ "${AREK_LEGACY_SOURCE_MODE:-}" != "isolated_research" ]]; then echo "Recovered historical entry; use canonical run_pipeline.py. To inspect native batching/diagnostics in a NEW output directory set AREK_LEGACY_SOURCE_MODE=isolated_research." >&2; exit 2; fi


###############################################################################
# SOLiD Color-Space Shotgun Metagenome Diagnostic & Preprocessing Preparation
# Project: PRJEB18041
# Data: ABI SOLiD 5500 single-end shotgun metagenome, gut microbiomes
#
# Output structure:
# PRJEB18041_solid/
# ├── 00_manifest/
# ├── 01_colorspace_qc/
# ├── 02_tool_check/
# ├── 03_preprocess_plan/
# ├── logs/
# └── status/
#
# NOTE: This script performs DIAGNOSTIC ONLY. No taxonomy profiling is done.
# Color-space reads require special handling:
# - Cannot use standard base-space tools (Kraken2, Bracken, bowtie2)
# - Requires color-space-aware mappers (SAET, BFAST, SHRiMP, bowtie with -C)
# - May need conversion to base-space via error correction (SAET) first
#
###############################################################################

############################
# User parameters
############################

PROJECT="PRJEB18041"

RAW_ROOT="${RAW_DATA_ROOT}/${PROJECT}/fastq"
OUT_ROOT="${RESULT_ROOT}/${PROJECT}_solid"

THREADS="${THREADS:-24}"

############################
# Output directories
############################

MANIFEST_DIR="${OUT_ROOT}/00_manifest"
QC_DIR="${OUT_ROOT}/01_colorspace_qc"
TOOL_CHECK_DIR="${OUT_ROOT}/02_tool_check"
PREPROCESS_PLAN_DIR="${OUT_ROOT}/03_preprocess_plan"
LOG_DIR="${OUT_ROOT}/logs"
STATUS_DIR="${OUT_ROOT}/status"

mkdir -p \
  "${MANIFEST_DIR}" \
  "${QC_DIR}" \
  "${TOOL_CHECK_DIR}" \
  "${PREPROCESS_PLAN_DIR}" \
  "${LOG_DIR}" \
  "${STATUS_DIR}"

MANIFEST="${MANIFEST_DIR}/manifest.tsv"
QC_REPORT="${QC_DIR}/colorspace_qc.tsv"
TOOL_CHECK="${TOOL_CHECK_DIR}/tool_check.tsv"
BOWTIE_COLORSPACE_CHECK="${TOOL_CHECK_DIR}/bowtie_colorspace_check.txt"
STATUS_FILE="${STATUS_DIR}/PRJEB18041_solid_status.txt"
LOG_FILE="${LOG_DIR}/prepare_PRJEB18041_solid.log"

############################
# Logging
############################

exec 1> >(tee -a "${LOG_FILE}")
exec 2>&1

log_info() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $*"
}

log_warn() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARNING] $*"
}

log_error() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2
}

############################
# Checks
############################

log_info "=================================================================="
log_info "SOLiD Color-Space Diagnostic Script"
log_info "Project: ${PROJECT}"
log_info "RAW_ROOT: ${RAW_ROOT}"
log_info "OUT_ROOT: ${OUT_ROOT}"
log_info "THREADS: ${THREADS}"
log_info "=================================================================="

if [[ ! -d "${RAW_ROOT}" ]]; then
  log_error "RAW_ROOT does not exist: ${RAW_ROOT}"
  exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
  log_error "python3 not found"
  exit 1
fi

###############################################################################
# 00_manifest: scan FASTQ files
###############################################################################

log_info "Scanning FASTQ files in ${RAW_ROOT}..."

python3 << PYSCN > "${MANIFEST}"
from pathlib import Path
import os

raw_root = Path("${RAW_ROOT}")

fastq_exts = [".fastq.gz", ".fq.gz", ".fastq", ".fq"]

def strip_fastq_ext(name):
    lower = name.lower()
    for ext in fastq_exts:
        if lower.endswith(ext):
            return name[:-len(ext)]
    return name

# Write header
print("sample_id\tfastq_path\tfile_size\tlayout\tread_space_guess")

sample_files = {}

for p in sorted(raw_root.rglob("*")):
    if not p.is_file():
        continue
    
    name_lower = p.name.lower()
    
    if p.name.endswith(".part"):
        continue
    
    if not any(name_lower.endswith(ext) for ext in fastq_exts):
        continue
    
    if p.stat().st_size == 0:
        continue
    
    sample_id = strip_fastq_ext(p.name)
    size = p.stat().st_size
    
    sample_files.setdefault(sample_id, []).append((str(p), size))

for sample, files in sorted(sample_files.items()):
    if len(files) == 1:
        path, size = files[0]
        print(f"{sample}\t{path}\t{size}\tSE\tunknown")
    else:
        for path, size in files:
            print(f"{sample}\t{path}\t{size}\tSE_multi\tunknown")
PYSCN

log_info "Manifest written: ${MANIFEST}"

###############################################################################
# 01_colorspace_qc: detect color-space characteristics
###############################################################################

log_info "Detecting color-space characteristics..."

python3 << PYQC > "${QC_REPORT}"
from pathlib import Path
import gzip
import sys
import re

def open_fastq(fpath):
    if fpath.endswith('.gz'):
        return gzip.open(fpath, 'rt', errors='replace')
    else:
        return open(fpath, 'r', errors='replace')

def read_fastq_record(f):
    header = f.readline().rstrip('\n')
    seq = f.readline().rstrip('\n')
    plus = f.readline().rstrip('\n')
    qual = f.readline().rstrip('\n')
    if not header:
        return None
    return (header, seq, plus, qual)

def is_base_space_read(seq):
    s = seq.upper()
    if not s:
        return False
    for c in s:
        if c not in 'ACGTN':
            return False
    return True

def is_color_space_read(seq):
    if not seq or len(seq) < 2:
        return False
    s = seq.upper()
    
    first = s[0]
    if first not in 'ACGTN':
        return False
    
    rest = s[1:]
    for c in rest:
        if c not in '01230.':
            return False
    
    return len(rest) > 0

def classify_read(seq):
    if is_base_space_read(seq):
        return 'base'
    elif is_color_space_read(seq):
        return 'color'
    else:
        return 'other'

manifest_file = Path("${MANIFEST}")

print("sample_id\tn_reads_checked\tread_len_min\tread_len_median\tread_len_max\t"
      "base_space_like\tcolor_space_like\tother_like\t"
      "dot_count_total\tdot_per_read_mean\tdot_read_fraction\t"
      "valid_color_fraction\texample_read_1\texample_read_2\texample_read_3")

if not manifest_file.exists():
    print("ERROR: Manifest not found", file=sys.stderr)
    sys.exit(1)

with open(manifest_file, 'r') as f:
    f.readline()
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 3:
            continue
        
        sample_id = parts[0]
        fastq_path = parts[1]
        
        if not Path(fastq_path).exists():
            continue
        
        try:
            fh = open_fastq(fastq_path)
        except Exception as e:
            print(f"ERROR: Cannot open {fastq_path}: {e}", file=sys.stderr)
            continue
        
        reads = []
        lengths = []
        base_count = 0
        color_count = 0
        other_count = 0
        dot_count_total = 0
        dot_read_count = 0
        examples = []
        
        for _ in range(10000):
            rec = read_fastq_record(fh)
            if not rec:
                break
            
            header, seq, plus, qual = rec
            if not seq:
                continue
            
            reads.append(seq)
            lengths.append(len(seq))
            
            classification = classify_read(seq)
            if classification == 'base':
                base_count += 1
            elif classification == 'color':
                color_count += 1
            else:
                other_count += 1
            
            dot_in_seq = seq.count('.') + seq.count('.')
            dot_count_total += dot_in_seq
            if '.' in seq or '.' in seq:
                dot_read_count += 1
            
            if len(examples) < 3:
                examples.append(seq[:80])
        
        fh.close()
        
        n_reads = len(reads)
        if n_reads == 0:
            continue
        
        lengths.sort()
        read_len_min = lengths[0] if lengths else 0
        read_len_max = lengths[-1] if lengths else 0
        read_len_median = lengths[n_reads // 2] if lengths else 0
        
        dot_per_read_mean = dot_count_total / n_reads if n_reads > 0 else 0.0
        dot_read_fraction = dot_read_count / n_reads if n_reads > 0 else 0.0
        
        # valid_color_fraction: reads that are classified as color-space
        valid_color_fraction = color_count / n_reads if n_reads > 0 else 0.0
        
        # Pad examples to 3
        while len(examples) < 3:
            examples.append("")
        
        output_parts = [
            sample_id,
            str(n_reads),
            str(read_len_min),
            str(read_len_median),
            str(read_len_max),
            str(base_count),
            str(color_count),
            str(other_count),
            str(dot_count_total),
            f"{dot_per_read_mean:.4f}",
            f"{dot_read_fraction:.4f}",
            f"{valid_color_fraction:.4f}",
            examples[0],
            examples[1],
            examples[2],
        ]
        
        print('\t'.join(output_parts))

PYQC

log_info "QC report written: ${QC_REPORT}"

###############################################################################
# 02_tool_check: check for color-space-aware tools
###############################################################################

log_info "Checking for color-space-aware tools..."

python3 << PYTOOL > "${TOOL_CHECK}"
import subprocess
from pathlib import Path
import shutil

tools_to_check = [
    'SAET',
    'BFAST',
    'bfast',
    'SHRiMP',
    'gmapper',
    'bowtie',
    'bowtie-build',
    'cutadapt',
    'metaphlan',
    'kraken2',
    'bracken',
    'bowtie2',
    'python',
    'perl',
]

print("tool\tfound\tpath\tversion_or_help_first_line")

for tool in tools_to_check:
    path = shutil.which(tool)
    if path:
        try:
            result = subprocess.run(
                [tool, '--version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            version_line = result.stdout.strip().split('\n')[0] if result.stdout else ''
            if not version_line:
                result = subprocess.run(
                    [tool, '--help'],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                version_line = result.stdout.strip().split('\n')[0] if result.stdout else ''
            version_line = version_line.replace('\t', ' ')[:100]
        except Exception as e:
            version_line = str(e)[:100]
        
        print(f"{tool}\tYES\t{path}\t{version_line}")
    else:
        print(f"{tool}\tNO\t\t")

PYTOOL

log_info "Tool check written: ${TOOL_CHECK}"

###############################################################################
# Check bowtie for color-space support
###############################################################################

log_info "Checking bowtie for color-space support..."

{
  echo "=== Checking bowtie --help for color-space mentions ==="
  if command -v bowtie >/dev/null 2>&1; then
    if bowtie --help 2>&1 | grep -i "color" ; then
      echo "Found color-space mention in bowtie --help"
    else
      echo "No color-space mention in bowtie --help"
    fi
  else
    echo "bowtie not found in PATH"
  fi
} > "${BOWTIE_COLORSPACE_CHECK}" 2>&1

{
  echo ""
  echo "=== Checking bowtie-build --help for color-space mentions ==="
  if command -v bowtie-build >/dev/null 2>&1; then
    if bowtie-build --help 2>&1 | grep -i "color" ; then
      echo "Found color-space mention in bowtie-build --help"
    else
      echo "No color-space mention in bowtie-build --help"
    fi
  else
    echo "bowtie-build not found in PATH"
  fi
} >> "${BOWTIE_COLORSPACE_CHECK}" 2>&1

log_info "Bowtie color-space check written: ${BOWTIE_COLORSPACE_CHECK}"

###############################################################################
# Generate status file
###############################################################################

log_info "Generating status file..."

cat > "${STATUS_FILE}" <<'EOFSTATUS'
================================================================================
PRJEB18041 - SOLiD Color-Space Shotgun Metagenome
Diagnostic & Preprocessing Preparation Status
================================================================================

PROJECT OVERVIEW
================================================================================
Project:        PRJEB18041
Data Type:      ABI SOLiD 5500 shotgun metagenome (single-end, color-space)
Source:         Gut microbiomes of patients with alcoholic dependence syndrome 
                and alcoholic liver cirrhosis
Data Location:  ${RAW_DATA_ROOT}/PRJEB18041/fastq
Output Location: ${RESULT_ROOT}/PRJEB18041_solid/

IMPORTANT: COLOR-SPACE READS
================================================================================
This dataset contains ABI SOLiD color-space reads, NOT base-space (A/C/G/T/N).

Color-space format:
  - First position: base (A/C/G/T/N)
  - Remaining positions: color values (0/1/2/3/.)
  - Example: T010..2.0.....21.1.31023...

This is fundamentally different from standard Illumina base-space reads.

NOT COMPATIBLE WITH
================================================================================
The following standard tools CANNOT be used directly on color-space reads:
  ✗ QIIME2 / VSEARCH (16S rRNA analysis)
  ✗ fastp (standard QC tool)
  ✗ bowtie2 (base-space mapper)
  ✗ Kraken2 (base-space taxonomic profiler)
  ✗ Bracken (base-space abundance estimator)

Attempting to process color-space reads as base-space will produce:
  - Incorrect quality control metrics
  - Meaningless taxonomic assignments
  - Completely unreliable results

REQUIRED FOR PROCESSING
================================================================================
To process PRJEB18041 color-space reads, you need ONE of:

1. Color-Space Aware Mappers:
   - SAET (SOLiD Accuracy Evaluation Tool)
     Purpose: Error correction & base-space conversion
     Status: Check in 02_tool_check/tool_check.tsv
   
   - BFAST (BLAT-like Fast Accurate Search Tool)
     Purpose: Direct color-space alignment
     Status: Check in 02_tool_check/tool_check.tsv
   
   - SHRiMP / gmapper
     Purpose: Direct color-space alignment
     Status: Check in 02_tool_check/tool_check.tsv
   
   - Original bowtie (with -C flag)
     Purpose: Can handle color-space with proper index
     Note: NOT bowtie2 (which is base-space only)
     Status: Check in 02_tool_check/bowtie_colorspace_check.txt

2. Two-Step Workflow (Recommended):
   Step 1: Use SAET to convert color-space → base-space (error correction)
   Step 2: Use standard base-space tools (fastp, Kraken2, Bracken, etc.)

DIAGNOSTIC OUTPUT
================================================================================
This script has generated:

1. 00_manifest/manifest.tsv
   - Sample identifiers and FASTQ file paths
   - Detected read layout (SE)
   
2. 01_colorspace_qc/colorspace_qc.tsv
   - Color-space vs. base-space classification
   - Read length statistics
   - Dot notation frequency (marks uncertain positions)
   - Example reads for manual inspection
   
3. 02_tool_check/tool_check.tsv
   - Installed tools on this system
   - Paths and versions
   - Check for SAET, BFAST, SHRiMP, bowtie, etc.
   
4. 02_tool_check/bowtie_colorspace_check.txt
   - Specific check for bowtie color-space support (-C flag)

NEXT STEPS
================================================================================
Option A: If SAET is available
  1. Run SAET correction: color-space → base-space with quality scores
  2. Re-process corrected reads as base-space (standard Kraken2/Bracken pipeline)
  3. Generate taxonomy tables
  
Option B: If BFAST/SHRiMP is available
  1. Build color-space index for reference database
  2. Align reads using color-space mapper
  3. Extract alignment-based taxonomy
  
Option C: If only base-space tools available
  1. STOP - Cannot proceed reliably
  2. Acquire SAET or contact data provider for base-space version
  3. Document as "not processable with current tools"

WARNING
================================================================================
DO NOT proceed with standard base-space tools without:
  ✗ Color-space → base-space conversion (SAET)
  ✗ OR color-space-aware alignment (BFAST/SHRiMP/bowtie -C)

Results will be meaningless and misleading.

For more information on SOLiD color-space format, refer to:
  - Applied Biosystems SOLiD Technical Notes
  - SAET documentation: https://www.ncbi.nlm.nih.gov/pubmed/
  - BFAST paper: https://www.ncbi.nlm.nih.gov/pubmed/

================================================================================
Generated: $(date)
================================================================================
EOFSTATUS

log_info "Status file written: ${STATUS_FILE}"

###############################################################################
# Summary
###############################################################################

log_info "=================================================================="
log_info "PRJEB18041 Diagnostic Complete"
log_info ""
log_info "Output files generated:"
log_info "  Manifest:    ${MANIFEST}"
log_info "  QC Report:   ${QC_REPORT}"
log_info "  Tool Check:  ${TOOL_CHECK}"
log_info "  Bowtie Info: ${BOWTIE_COLORSPACE_CHECK}"
log_info "  Status:      ${STATUS_FILE}"
log_info ""
log_info "IMPORTANT:"
log_info "  - This project contains SOLiD COLOR-SPACE reads (not base-space)"
log_info "  - Cannot use standard fastp/Kraken2/Bracken/bowtie2"
log_info "  - Requires SAET, BFAST, SHRiMP, or bowtie with -C flag"
log_info "  - Review status file for next steps"
log_info "=================================================================="

log_info "Diagnostic complete. Log written: ${LOG_FILE}"
