# =========================================================
# PCA of RNA-seq data — WW and WD scenarios
# Produces 3 plots for WW and 5 plots for WD:
#   - WW / WD coloured by genotype
#   - WW / WD coloured by SWP mean drought (7bar, 5bar)
#   - WD only: coloured by days under drought threshold (7bar, 5bar)
# Note: days_under_threshold plots are WD-only — WW plants
#       were never under water deficit so those values are absent
# =========================================================

library(tximport)
library(readr)
library(DESeq2)
library(readxl)
library(GenomicFeatures)
library(ggplot2)
library(viridis)
library(dplyr)
library(AnnotationDbi)

rm(list = ls())

# =========================================================
# Shared setup
# =========================================================
samples      <- read_excel("sample_metadata.xlsx")
kallisto_dir <- "kallisto_results"
gtf_file     <- "GCF_036512215.1_SLM_r2.1_genomic.gtf"
dir.create("grafici_PCA", showWarnings = FALSE)

# Build transcript -> gene mapping (done once, shared across all PCAs)
txdb    <- makeTxDbFromGFF(gtf_file, format = "gtf")
tx2gene <- AnnotationDbi::select(txdb,
                                 keys    = keys(txdb, keytype = "TXNAME"),
                                 keytype = "TXNAME",
                                 columns = "GENEID")
colnames(tx2gene) <- c("TXNAME", "GENEID")

# =========================================================
# Function: build VST-transformed DESeq2 object for a scenario
# =========================================================
build_vsd <- function(samples_subset) {
  files <- file.path(kallisto_dir, samples_subset$sample, "abundance.tsv")
  names(files) <- samples_subset$sample
  stopifnot(all(file.exists(files)))

  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
  dds <- DESeqDataSetFromTximport(txi, colData = samples_subset, design = ~ geno)

  # Pre-filter: keep genes with >=5 reads in >=2 samples
  keep <- rowSums(counts(dds) >= 5) >= 2
  dds  <- dds[keep, ]

  # blind = TRUE: transformation unaware of design, appropriate for QC/exploration
  vst(dds, blind = TRUE)
}

# =========================================================
# Function: compute PCA and optionally flip PC1
# (WT 0_0_0_0 kept on positive side for visual consistency)
# =========================================================
compute_pca <- function(vsd) {
  pca_data   <- plotPCA(vsd, intgroup = "geno", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

  wt_mean <- mean(pca_data$PC1[pca_data$geno == "0_0_0_0"])
  if (wt_mean < 0) pca_data$PC1 <- -pca_data$PC1

  list(pca_data = pca_data, percentVar = percentVar)
}

# =========================================================
# Function: plot PCA coloured by genotype (discrete)
# =========================================================
plot_pca_geno <- function(pca_data, percentVar, scenario) {
  ggplot(pca_data, aes(PC1, PC2, color = geno)) +
    geom_point(size = 3) +
    labs(
      x     = paste0("PC1: ", percentVar[1], "% variance"),
      y     = paste0("PC2: ", percentVar[2], "% variance"),
      title = paste0("PCA of RNA-seq samples (", scenario, " scenario)")
    ) +
    scale_color_manual(values = c(
      `0_0_0_0`       = "#000000",
      `42_5_5|2_5`    = "#BF6193",
      `42|35_5_5|2_5` = "#FDAE8A",
      `35_5_5|2_5`    = "#6A4D8D"
    )) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(hjust = 0.5)
    )
}

# =========================================================
# Function: plot PCA coloured by a continuous variable
# Used for both SWP and days_under_threshold
# =========================================================
plot_pca_continuous <- function(pca_data, percentVar, color_var, scenario,
                                legend_label, viridis_option, plot_title) {
  ggplot(pca_data, aes(PC1, PC2, color = .data[[color_var]], shape = geno)) +
    geom_point(size = 4) +
    scale_color_viridis_c(option = viridis_option, name = legend_label) +
    labs(
      x     = paste0("PC1: ", percentVar[1], "% variance"),
      y     = paste0("PC2: ", percentVar[2], "% variance"),
      title = plot_title
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(hjust = 0.5)
    )
}

# =========================================================
# Run for both scenarios
# =========================================================
for (scenario in c("WW", "WD")) {

  samples_sub <- samples %>%
    filter(scenario == !!scenario) %>%
    as.data.frame()
  rownames(samples_sub) <- samples_sub$sample

  # Ensure numeric columns
  samples_sub$SWP_mean_Drought_7bar             <- as.numeric(samples_sub$SWP_mean_Drought_7bar)
  samples_sub$SWP_mean_Drought_5bar             <- as.numeric(samples_sub$SWP_mean_Drought_5bar)
  samples_sub$days_under_threshold_Drought_7bar <- as.numeric(samples_sub$days_under_threshold_Drought_7bar)
  samples_sub$days_under_threshold_Drought_5bar <- as.numeric(samples_sub$days_under_threshold_Drought_5bar)
  samples_sub$geno <- factor(samples_sub$geno)

  vsd        <- build_vsd(samples_sub)
  pca        <- compute_pca(vsd)
  pca_data   <- pca$pca_data
  percentVar <- pca$percentVar

  # Merge all continuous metadata back into pca_data for colouring
  pca_data <- pca_data %>%
    left_join(
      samples_sub %>% select(sample,
                             SWP_mean_Drought_7bar,
                             SWP_mean_Drought_5bar,
                             days_under_threshold_Drought_7bar,
                             days_under_threshold_Drought_5bar),
      by = c("name" = "sample")
    )

  # --- Plot 1: coloured by genotype ---
  p_geno <- plot_pca_geno(pca_data, percentVar, scenario)
  ggsave(paste0("grafici_PCA/PCA_", scenario, "_geno.svg"), p_geno, width = 8, height = 6)

  # --- Plot 2: coloured by SWP 7bar ---
  p_swp7 <- plot_pca_continuous(pca_data, percentVar,
                                color_var      = "SWP_mean_Drought_7bar",
                                scenario       = scenario,
                                legend_label   = "SWP < 7bar",
                                viridis_option = "plasma",
                                plot_title     = paste0(scenario, " — PCA coloured by SWP mean drought 7bar"))
  ggsave(paste0("grafici_PCA/PCA_", scenario, "_SWP_7bar.svg"), p_swp7, width = 8, height = 6)

  # --- Plot 3: coloured by SWP 5bar ---
  p_swp5 <- plot_pca_continuous(pca_data, percentVar,
                                color_var      = "SWP_mean_Drought_5bar",
                                scenario       = scenario,
                                legend_label   = "SWP < 5bar",
                                viridis_option = "magma",
                                plot_title     = paste0(scenario, " — PCA coloured by SWP mean drought 5bar"))
  ggsave(paste0("grafici_PCA/PCA_", scenario, "_SWP_5bar.svg"), p_swp5, width = 8, height = 6)

  # --- Plots 4-5: days under threshold — WD only ---
  # WW plants were never under water deficit, so this variable is absent
  if (scenario == "WD") {

    p_days7 <- plot_pca_continuous(pca_data, percentVar,
                                   color_var      = "days_under_threshold_Drought_7bar",
                                   scenario       = scenario,
                                   legend_label   = "Days < 7bar",
                                   viridis_option = "plasma",
                                   plot_title     = "WD — PCA coloured by days under drought threshold (7bar)")
    ggsave("grafici_PCA/PCA_WD_days_7bar.svg", p_days7, width = 8, height = 6)

    p_days5 <- plot_pca_continuous(pca_data, percentVar,
                                   color_var      = "days_under_threshold_Drought_5bar",
                                   scenario       = scenario,
                                   legend_label   = "Days < 5bar",
                                   viridis_option = "magma",
                                   plot_title     = "WD — PCA coloured by days under drought threshold (5bar)")
    ggsave("grafici_PCA/PCA_WD_days_5bar.svg", p_days5, width = 8, height = 6)
  }

  message("PCA plots saved for scenario: ", scenario)
}
