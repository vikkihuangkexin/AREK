#!/usr/bin/env bash
set -euo pipefail
# Exact biological filtering flags from run_PRJNA613834_metagenome_taxonomy_v2.sh.
index=$1;threads=$2;sort_memory=$3;r1=$4;r2=$5;bam=$6;o1=$7;o2=$8
if [[ -n "$r2" ]]; then
bowtie2 -x "$index" -1 "$r1" -2 "$r2" --very-sensitive --threads "$threads" |
 samtools view -@ "$threads" -b -f 12 -F 256 - |
 samtools sort -@ "$threads" -m "$sort_memory" -n -o "$bam" -
samtools quickcheck -v "$bam"
samtools fastq -@ "$threads" -1 "$o1" -2 "$o2" -0 /dev/null -s /dev/null -n "$bam"
else
 bowtie2 -x "$index" -U "$r1" --very-sensitive --threads "$threads" |
  samtools view -@ "$threads" -b -f 4 -F 256 - > "$bam"
 samtools quickcheck -v "$bam"
 samtools fastq -@ "$threads" "$bam" > "$o1"
fi
