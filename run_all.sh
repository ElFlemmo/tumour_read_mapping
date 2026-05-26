#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 <fastq_dir> <reference.fasta> <annotation.gff3>"
    exit 1
}

[[ $# -lt 3 ]] && usage

FASTQ_DIR=$1
FASTA=$2
GFF=$3
INDEX_PREFIX=$(basename "${FASTA%.*}")
SCRIPT_DIR=$(dirname "$0")

echo "Building Bowtie2 index from $FASTA"
bowtie2-build "$FASTA" "$INDEX_PREFIX"

for R1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    bash "${SCRIPT_DIR}/read_mapping.sh" "$SAMPLE" "$FASTQ_DIR" "$INDEX_PREFIX" "$GFF"
done

echo "All samples processed."
