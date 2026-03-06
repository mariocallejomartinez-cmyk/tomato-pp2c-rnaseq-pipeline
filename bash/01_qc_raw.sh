#!/bin/bash
# ==============================================================
# QC on raw reads using FastQC and MultiQC
# Run from rna_seq/ directory
# Input: samples.txt (one sample name per line)
# Output: QC_1/new/<sample>/
# ==============================================================

FASTQC="../../bin/FastQC/fastqc"
SAMPLE_FILE="samples.txt"

while read sample; do
    mkdir -p QC_1/new/${sample}
    ${FASTQC} raw_data/${sample}_1.fq.gz -o QC_1/new/${sample}/ -t 32
    ${FASTQC} raw_data/${sample}_2.fq.gz -o QC_1/new/${sample}/ -t 32
done < "${SAMPLE_FILE}"

multiqc QC_1/new/ -o QC_1/new/
