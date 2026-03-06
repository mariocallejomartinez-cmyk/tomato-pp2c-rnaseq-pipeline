#!/bin/bash
# ==============================================================
# Kallisto quantification on trimmed paired-end reads
# Run from project root directory
# Input: samples.txt (one sample name per line)
#        trimmed reads from fastp (02_trimming_qc.sh output)
# Output: rna_seq/align_quant/kallisto_results/<sample>/
# Note: bootstraps set to 0 — not required for DESeq2
# ==============================================================

INDEX=rna_seq/align_quant/transcripts.idx
READS_DIR=rna_seq/trimmed/
OUTPUT_DIR=rna_seq/align_quant/kallisto_results
KALLISTO=../bin/kallisto/kallisto
SAMPLE_FILE="samples.txt"

mkdir -p ${OUTPUT_DIR}

while read sample; do
    R1=${READS_DIR}/${sample}_trimmed/${sample}_trimmed_1.fq.gz
    R2=${READS_DIR}/${sample}_trimmed/${sample}_trimmed_2.fq.gz

    ${KALLISTO} quant \
        -i ${INDEX} \
        -o ${OUTPUT_DIR}/${sample} \
        -b 0 \
        -t 32 \
        ${R1} ${R2}

done < "${SAMPLE_FILE}"
