#cheese
#!/usr/bin/env bash
set -euo pipefail

SAMPLE=$1
FASTQ_DIR=$2
INDEX=$3
GFF=$4

echo "Running fastp for $SAMPLE"
fastp \
  -i "${FASTQ_DIR}/${SAMPLE}_1.fastq" \
  -I "${FASTQ_DIR}/${SAMPLE}_2.fastq" \
  -o "${SAMPLE}_1.cleaned.fastq" \
  -O "${SAMPLE}_2.cleaned.fastq" \
  --detect_adapter_for_pe \
  --thread 48 \
  --html "${SAMPLE}_fastp.html" \
  --json "${SAMPLE}_fastp.json"

echo "Aligning $SAMPLE with Bowtie2"
bowtie2 -p 48 -x "$INDEX" \
  -1 "${SAMPLE}_1.cleaned.fastq" \
  -2 "${SAMPLE}_2.cleaned.fastq" | \
  samtools view -@ 48 -bS -F 4 - | samtools sort -@ 48 -o "${SAMPLE}.sorted.bam"

samtools index "${SAMPLE}.sorted.bam"

echo "Counting reads for $SAMPLE"
featureCounts -a "$GFF" -o "${SAMPLE}.cds.counts.txt" \
  -t CDS -g ID -M -O -s 0 -p -T 48 "${SAMPLE}.sorted.bam"

echo "Done with $SAMPLE"
