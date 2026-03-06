# Functional and Molecular Characterization of Tomato PP2C Mutants — RNA-seq Analysis Pipeline

This repository contains the bioinformatics pipelines developed for my BSc thesis in Biotechnology (University of Milan, 2025). The study investigates the transcriptomic response to drought stress in *Solanum lycopersicum* PP2C mutants using RNA-seq data.

## Biological context

Four PP2C mutant genotypes with deletions in four genes and a wild-type control were grown under two irrigation scenarios: well-watered (WW) and water deficit (WD). RNA-seq was performed on all samples and analysed to characterize genotype-dependent and scenario-dependent gene expression changes.

## Repository structure

```
thesis-pipelines/
├── README.md
├── bash/
│   ├── 01_qc_raw.sh           # Quality control on raw reads (FastQC + MultiQC)
│   ├── 02_trimming_qc.sh      # Trimming with fastp + QC on trimmed reads (FastQC + MultiQC)
│   ├── 03_kallisto_quant.sh   # Pseudoalignment and quantification with Kallisto
│   └── samples.txt            # One sample name per line — edit before running
└── R/
    ├── 00_build_annotation_db.R    # Build org.Slycopersicum.eg.db from NCBI (run once)
    ├── 01_pca.R                    # Exploratory PCA (WW and WD, multiple colourings)
    ├── 02_deseq2.R                 # Differential expression analysis by genotype
    ├── 03_venn.R                   # Venn diagrams of DEGs across genotype contrasts
    ├── 04_gprofiler_rrvgo.R        # GO/KEGG enrichment + rrvgo redundancy reduction
    └── 05_target_genes_scenario.R  # Scenario effect on target PP2C genes in WT
```

## Workflow overview

```
Raw reads
    │
    ▼
01_qc_raw.sh          → FastQC + MultiQC report on raw reads
    │
    ▼
02_trimming_qc.sh     → fastp trimming + FastQC + MultiQC on trimmed reads
    │
    ▼
03_kallisto_quant.sh  → Kallisto quantification per sample
    │
    ▼
01_pca.R              → Exploratory PCA (WW and WD scenarios)
    │
    ▼
02_deseq2.R           → DEG analysis by genotype (WW and WD separately)
    │
    ├──▶ 03_venn.R                  → Venn diagrams of DEGs across contrasts
    │
    └──▶ 04_gprofiler_rrvgo.R       → Functional enrichment on common DEGs (WW)
         05_target_genes_scenario.R → Scenario effect on PP2C target genes (WT only)
```

## Requirements

### Bash

All bash scripts were run on a Linux server. The following tools are required and are expected at the relative paths defined in each script (adjust if needed):

- FastQC
- MultiQC
- fastp
- Kallisto (bootstraps set to 0 — not required for DESeq2)

Tool versions were not recorded.

### R

R version 4.5.1 (2025-06-13), running on macOS Sequoia 15.6 (aarch64).

Key packages:

| Package | Version |
|---|---|
| DESeq2 | 1.48.2 |
| tximport | 1.36.1 |
| GenomicFeatures | 1.60.0 |
| AnnotationDbi | 1.70.0 |
| gprofiler2 | 0.2.4 |
| rrvgo | 1.20.0 |
| org.Slycopersicum.eg.db | 1.0 |
| ggplot2 | 4.0.1 |
| VennDiagram | 1.7.3 |
| apeglm | 1.30.0 |
| openxlsx | 4.2.8.1 |
| viridis | 0.6.5 |

The `org.Slycopersicum.eg.db` annotation package must be built locally before running enrichment analysis. See `00_build_annotation_db.R`. Note: `makeOrgPackageFromNCBI()` may fail on macOS due to URL handling issues — run the script on Windows if this occurs.

## Usage

### Bash pipelines

1. Populate `bash/samples.txt` with one sample name per line.
2. Adjust the paths to tool binaries at the top of each script if needed.
3. Run scripts sequentially from the `rna_seq/` directory:

```bash
bash 01_qc_raw.sh
bash 02_trimming_qc.sh
bash 03_kallisto_quant.sh
```

### R scripts

1. Run `00_build_annotation_db.R` once to install the tomato annotation package.
2. Place `sample_metadata.xlsx` and the Kallisto results folder in the working directory.
3. Run scripts in order (01 → 05).

`sample_metadata.xlsx` must contain at minimum the following columns: `sample`, `geno`, `scenario`, `SWP_mean_Drought_7bar`, `SWP_mean_Drought_5bar`, `days_under_threshold_Drought_7bar`, `days_under_threshold_Drought_5bar`.

## Key analysis notes

- **Pre-filtering:** PCA uses ≥5 reads in ≥2 samples; DEG analysis and enrichment background both use ≥10 reads in ≥4 samples to avoid p-value distribution artefacts caused by low-count genes.
- **DEG thresholds:** padj < 0.05, |log2FC| ≥ 1 (≥2-fold change), with apeglm log2FC shrinkage.
- **Enrichment:** run on common DEGs across all genotype contrasts without up/down separation, due to low gene count. Background defined as all genes passing the DESeq2 pre-filter, mapped to Entrez IDs via the NCBI GTF.
- **Genome annotation:** GCF_036512215.1_SLM_r2.1 (*Solanum lycopersicum*).
- **PC1 orientation:** flipped where necessary so that WT (0_0_0_0) projects consistently on the positive side across all PCA plots.

## Author

Mario Callejo Martinez — BSc Biotechnology, University of Milan, 2025
