#cheese
#!/usr/bin/env bash
set -euo pipefail

SAMPLE=$1
FASTQ_DIR=$2
INDEX=$3
GFF=$4

echo "Running fastp for $SAMPLE"
fastp \
  -i "${FASTQ_DIR}/${SAMPLE}_1.fastq.gz" \
  -I "${FASTQ_DIR}/${SAMPLE}_2.fastq.gz" \
  -o "${SAMPLE}_1.cleaned.fastq.gz" \
  -O "${SAMPLE}_2.cleaned.fastq.gz" \
  --detect_adapter_for_pe \
  --thread 16 \
  --html "${SAMPLE}_fastp.html" \
  --json "${SAMPLE}_fastp.json"

echo "Aligning $SAMPLE with Bowtie2"
bowtie2 -x "$INDEX" \
  -1 "${SAMPLE}_1.cleaned.fastq.gz" \
  -2 "${SAMPLE}_2.cleaned.fastq.gz" | \
  samtools view -bS -F 4 - | samtools sort -o "${SAMPLE}.sorted.bam"

samtools index "${SAMPLE}.sorted.bam"

echo "Counting reads for $SAMPLE"
featureCounts -a "$GFF" -o "${SAMPLE}.cds.counts.txt" \
  -t CDS -g ID -M -O -s 0 -p "${SAMPLE}.sorted.bam"

echo "Done with $SAMPLE"
