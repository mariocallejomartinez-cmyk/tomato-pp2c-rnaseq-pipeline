# ============================================================
# Scenario effect (WD vs WW) on target PP2C genes in WT
# Design: ~ scenario (WT samples only)
# Reference: WW (control irrigation)
# Target genes: PP2C-2, LOC101248778, ABI2, PP2C-1
# ============================================================

rm(list = ls())

library(readr)
library(dplyr)
library(readxl)
library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(openxlsx)
library(ggplot2)
library(viridis)
library(AnnotationDbi)

# ============================================================
# 1) Sample metadata — WT only
# ============================================================
wt_geno      <- "0_0_0_0"
ref_scenario <- "WW"

samples <- read_excel("sample_metadata.xlsx")

WT_samples <- samples %>%
  filter(geno == wt_geno) %>%
  as.data.frame()
rownames(WT_samples) <- WT_samples$sample

# ============================================================
# 2) Kallisto abundance files
# ============================================================
kallisto_dir <- "kallisto_results"
files <- file.path(kallisto_dir, WT_samples$sample, "abundance.tsv")
names(files) <- WT_samples$sample
stopifnot(all(file.exists(files)))

# ============================================================
# 3) Transcript -> gene mapping from GTF
# ============================================================
gtf_file <- "GCF_036512215.1_SLM_r2.1_genomic.gtf"
txdb     <- makeTxDbFromGFF(gtf_file, format = "gtf")
tx2gene  <- AnnotationDbi::select(txdb,
                                  keys    = keys(txdb, keytype = "TXNAME"),
                                  keytype = "TXNAME",
                                  columns = "GENEID")
colnames(tx2gene) <- c("TXNAME", "GENEID")

# ============================================================
# 4) Import kallisto output and build DESeq2 object
# ============================================================
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)

dds <- DESeqDataSetFromTximport(txi, colData = WT_samples, design = ~ scenario)

# Pre-filter: keep genes with >=5 reads in >=2 samples
keep <- rowSums(counts(dds) >= 5) >= 2
dds  <- dds[keep, ]

# Set WW as reference scenario
dds$scenario <- relevel(dds$scenario, ref = ref_scenario)

dds <- DESeq(dds)
resultsNames(dds)

# ============================================================
# 5) Extract results for target PP2C genes
# ============================================================
target_genes <- c("PP2C-2", "LOC101248778", "ABI2", "PP2C-1")

res_wt <- results(dds, contrast = c("scenario", "WD", "WW"))

res_target <- as.data.frame(res_wt[target_genes, ])
res_target$gene <- rownames(res_target)

dir.create("DEG_results/effetto_scenario", showWarnings = FALSE, recursive = TRUE)
write.xlsx(res_target,
           file = "DEG_results/effetto_scenario/target_genes_WDvsWW_WT.xlsx",
           rowNames = FALSE)

# ============================================================
# 6) Barplot — log2FC per target gene (WD vs WW in WT)
# ============================================================
my_colors <- viridis(nrow(res_target), option = "magma", begin = 0.2, end = 0.9)
res_target$color       <- my_colors
res_target$significant <- res_target$padj < 0.05

p <- ggplot(res_target, aes(x = gene, y = log2FoldChange,
                             fill = color, alpha = significant)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_identity() +
  labs(
    title = "Scenario effect on target PP2Cs in wild type (WD vs WW)",
    x     = "Target gene",
    y     = "Log2 fold change (WD / WW)",
    alpha = "Significant (padj < 0.05)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "bottom",
    plot.title      = element_text(hjust = 0.5)
  )

ggsave("DEG_results/effetto_scenario/plot_WDvsWW_inWT.png",
       p, width = 11, height = 8)
