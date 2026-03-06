# ============================================================
# Differential expression analysis — WW and WD scenarios
# Design: ~ geno (effect of genotype within each scenario)
# Reference genotype: 0_0_0_0 (WT)
# Output: per-contrast DEG tables + common DEGs across all contrasts
# Thresholds: padj < 0.05, |log2FC| >= 1 (>=2-fold change)
# ============================================================

library(tximport)
library(readr)
library(DESeq2)
library(readxl)
library(GenomicFeatures)
library(ggplot2)
library(writexl)
library(apeglm)
library(openxlsx)
library(dplyr)
library(AnnotationDbi)

rm(list = ls())

# ============================================================
# Shared setup
# ============================================================
samples      <- read_excel("sample_metadata.xlsx")
kallisto_dir <- "kallisto_results"
gtf_file     <- "GCF_036512215.1_SLM_r2.1_genomic.gtf"

padj_threshold   <- 0.05
log2fc_threshold <- 1  # |log2FC| >= 1 corresponds to >=2-fold change

# Build transcript -> gene mapping (done once, shared across both scenarios)
txdb    <- makeTxDbFromGFF(gtf_file, format = "gtf")
tx2gene <- AnnotationDbi::select(txdb,
                                 keys    = keys(txdb, keytype = "TXNAME"),
                                 keytype = "TXNAME",
                                 columns = "GENEID")
colnames(tx2gene) <- c("TXNAME", "GENEID")

# ============================================================
# Run DEG analysis for each scenario
# ============================================================
for (scenario in c("WW", "WD")) {

  message("--- Scenario: ", scenario, " ---")

  # Subset metadata
  samples_sub <- samples %>%
    filter(scenario == !!scenario) %>%
    as.data.frame()
  rownames(samples_sub) <- samples_sub$sample
  samples_sub$geno <- factor(samples_sub$geno)

  # Kallisto abundance files
  files <- file.path(kallisto_dir, samples_sub$sample, "abundance.tsv")
  names(files) <- samples_sub$sample
  stopifnot(all(file.exists(files)))

  # Import kallisto output
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
  stopifnot(all(colnames(txi$counts) == rownames(samples_sub)))

  # Build DESeq2 object
  dds <- DESeqDataSetFromTximport(txi, colData = samples_sub, design = ~ geno)
  dds$geno <- relevel(dds$geno, ref = "0_0_0_0")

 # Pre-filter: keep genes with >=10 reads in >=4 samples
 # (increased from default to avoid p-value distribution artefacts
 # caused by low-count genes, as observed in preliminary analysis)
  keep <- rowSums(counts(dds) >= 10) >= 4
  dds  <- dds[keep, ]

  # Run DESeq2
  dds <- DESeq(dds)
  message("Model coefficients: ", paste(resultsNames(dds), collapse = ", "))

  # Output directory
  outdir <- paste0("DEG_results_", scenario)
  dir.create(outdir, showWarnings = FALSE)

  # ----------------------------------------------------------
  # Extract and save DEGs — per contrast
  # ----------------------------------------------------------
  coef_names <- resultsNames(dds)

  for (coef_name in coef_names) {
    if (tolower(coef_name) == "intercept") next

    message("Processing: ", coef_name)

    res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
    res_df <- as.data.frame(res_shrunk)
    res_df$gene_id <- rownames(res_df)
    res_df <- res_df %>% select(gene_id, everything())

    sig_df <- res_df %>%
      filter(!is.na(padj)) %>%
      filter(padj < padj_threshold & abs(log2FoldChange) >= log2fc_threshold) %>%
      mutate(regulation = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
      arrange(desc(log2FoldChange))

    up_count   <- sum(sig_df$regulation == "Up",   na.rm = TRUE)
    down_count <- sum(sig_df$regulation == "Down", na.rm = TRUE)
    message("  Up: ", up_count, " | Down: ", down_count)

    if (nrow(sig_df) > 0) {
      out_file <- file.path(outdir, paste0("DEGs_", coef_name, ".xlsx"))
      write.xlsx(sig_df %>% select(gene_id, log2FoldChange, pvalue, padj, regulation),
                 file = out_file, rowNames = FALSE)
    } else {
      message("  No significant DEGs for: ", coef_name)
    }
  }

  # ----------------------------------------------------------
  # Identify DEGs common to all contrasts
  # ----------------------------------------------------------
  deg_files   <- list.files(outdir, pattern = "^DEGs_", full.names = TRUE)
  deg_lists   <- lapply(deg_files, function(f) read.xlsx(f)$gene_id)
  common_degs <- Reduce(intersect, deg_lists)

  message("Genes DE in all contrasts (", scenario, "): ", length(common_degs))

  if (length(common_degs) > 0) {
    write.xlsx(data.frame(gene_id = common_degs),
               file = file.path(outdir, "common_DEGs_all_contrasts.xlsx"),
               rowNames = FALSE)
  }
}
