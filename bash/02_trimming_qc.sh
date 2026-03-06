#!/bin/bash
# ==============================================================
# Trimming with fastp + QC on trimmed reads (FastQC + MultiQC)
# Run from rna_seq/ directory
# Input: samples.txt (one sample name per line)
# Output: trimmed/<sample>/, QC_2/new/<sample>/
# Parameters: q20, min length 40bp, PE adapter detection,
#             poly-G and poly-X trimming
# ==============================================================

FASTQC="../../bin/FastQC/fastqc"
FASTP="../../bin/fastp"
SAMPLE_FILE="samples.txt"

while read sample; do
    mkdir -p trimmed/${sample}_trimmed trimmed/json_reports trimmed/html_reports

    ${FASTP} \
        -i raw_data/${sample}_1.fq.gz \
        -I raw_data/${sample}_2.fq.gz \
        -q 20 -u 40 -l 40 \
        --detect_adapter_for_pe --trim_poly_g --trim_poly_x \
        --thread 16 \
        -o trimmed/${sample}_trimmed/${sample}_trimmed_1.fq.gz \
        -O trimmed/${sample}_trimmed/${sample}_trimmed_2.fq.gz \
        -j trimmed/json_reports/${sample}_report.json \
        -h trimmed/html_reports/${sample}_report.html

    mkdir -p QC_2/new/${sample}
    ${FASTQC} trimmed/${sample}_trimmed/${sample}_trimmed_1.fq.gz -o QC_2/new/${sample}/ -t 32
    ${FASTQC} trimmed/${sample}_trimmed/${sample}_trimmed_2.fq.gz -o QC_2/new/${sample}/ -t 32

done < "${SAMPLE_FILE}"

multiqc QC_2/new/ -o QC_2/new/
