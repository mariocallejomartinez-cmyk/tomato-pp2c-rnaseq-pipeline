# ============================================================
# Venn diagrams of DEGs by genotype effect
# WW: 3-way venn (3 mutant genotypes vs WT)
# WD: 2-way venn (2 mutant genotypes vs WT)
# Note: 35_5_5|2_5 absent in WD due to missing samples
# ============================================================

library(readxl)
library(VennDiagram)
library(viridis)

# ============================================================
# Load DEG tables
# ============================================================
deg1 <- read_excel("DEG_results/effetto_geno/WW/significant_DEGs_geno_35_5_5.2_5_vs_0_0_0_0.xlsx")
deg2 <- read_excel("DEG_results/effetto_geno/WW/significant_DEGs_geno_42_5_5.2_5_vs_0_0_0_0.xlsx")
deg3 <- read_excel("DEG_results/effetto_geno/WW/significant_DEGs_geno_42.35_5_5.2_5_vs_0_0_0_0.xlsx")
deg4 <- read_excel("DEG_results/effetto_geno/WD/significant_DEGs_geno_42.35_5_5.2_5_vs_0_0_0_0.xlsx")
deg5 <- read_excel("DEG_results/effetto_geno/WD/significant_DEGs_geno_42_5_5.2_5_vs_0_0_0_0.xlsx")

# ============================================================
# Gene lists
# ============================================================
gene_list_WW <- list(
  `35_5_5|2_5`    = deg1$GENEID,
  `42_5_5|2_5`    = deg2$GENEID,
  `42|35_5_5|2_5` = deg3$GENEID
)

gene_list_WD <- list(
  `42_5_5|2_5`    = deg5$GENEID,
  `42|35_5_5|2_5` = deg4$GENEID
)

dir.create("Venn/effetto_geno", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Venn WW — 3-way
# ============================================================
my_colors_WW <- viridis(length(gene_list_WW), option = "magma", begin = 0.2, end = 0.8)

venn.diagram(
  x        = gene_list_WW,
  filename = "Venn/effetto_geno/Common_deregulated_genes_WW.tiff",
  fill     = my_colors_WW,
  alpha    = 0.7,
  cex      = 1.4,
  cat.cex  = 1.2,
  cat.col  = my_colors_WW,
  res      = 600,
  main     = "Common deregulated genes WW"
)

# ============================================================
# Venn WD — 2-way
# ============================================================
my_colors_WD <- c("#BF6193", "#FDAE8A")

venn.diagram(
  x        = gene_list_WD,
  filename = "Venn/effetto_geno/Common_deregulated_genes_WD.tiff",
  fill     = my_colors_WD,
  alpha    = 0.7,
  cex      = 1.4,
  cat.cex  = 1,
  cat.col  = my_colors_WD,
  height   = 5,
  width    = 5,
  units    = "in",
  cat.pos  = c(330, 25),
  cat.dist = 0.05,
  res      = 600,
  main     = "Common deregulated genes WD"
)
