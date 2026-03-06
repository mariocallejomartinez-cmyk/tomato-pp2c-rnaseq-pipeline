# ============================================================
# Functional enrichment analysis with g:Profiler + rrvgo
# Input: common DEGs across all genotype contrasts (WW scenario)
# Background: all genes tested in DESeq2 (mapped to Entrez)
# Note: analysis run on common DEGs without up/down separation
#       due to low gene count — splitting produced 0-1 classes
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
library(stringr)
library(tidyr)
library(rtracklayer)
library(gprofiler2)
library(rrvgo)
library(GOSemSim)
library(org.Slycopersicum.eg.db)
library(AnnotationDbi)

rm(list = ls())

# ============================================================
# 1) Sample metadata
# ============================================================
samples <- read_excel("sample_metadata.xlsx")

WW_samples <- samples %>%
  filter(scenario == "WW") %>%
  as.data.frame()
rownames(WW_samples) <- WW_samples$sample

wt_geno      <- "0_0_0_0"
mutant_genos <- c("35_5_5|2_5", "42_5_5|2_5", "42|35_5_5|2_5")

# ============================================================
# 2) Kallisto abundance files
# ============================================================
kallisto_dir <- "kallisto_results"
files <- file.path(kallisto_dir, WW_samples$sample, "abundance.tsv")
names(files) <- WW_samples$sample
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

dds <- DESeqDataSetFromTximport(txi, colData = WW_samples, design = ~ geno)

# Pre-filter: >=10 reads in >=4 samples
# Note: more stringent than DEG analysis filter (>=5 in >=3)
# used here to define a high-confidence background for enrichment
keep <- rowSums(counts(dds) >= 10) >= 4
dds  <- dds[keep, ]

dds$geno <- relevel(dds$geno, ref = wt_geno)

# ============================================================
# 5) Build LOC -> Entrez mapping from GTF
# ============================================================
g <- import(gtf_file)
m <- as.data.frame(mcols(g))

gtf_genes <- if ("type" %in% colnames(m)) m[m$type == "gene", ] else m

mapping <- gtf_genes %>%
  as.data.frame() %>%
  mutate(db_xref = as.character(db_xref)) %>%
  mutate(Entrez  = str_extract(db_xref, "GeneID:\\d+")) %>%
  mutate(Entrez  = sub("GeneID:", "", Entrez)) %>%
  dplyr::select(gene_id, gene, gene_biotype, Entrez, everything()) %>%
  distinct()

write.csv(mapping, "NCBI_geneid_to_entrez_mapping.csv", row.names = FALSE)
message("Mapping table saved: LOC ID -> Entrez ID")

# ============================================================
# 6) Define background gene set
# ============================================================
background_genes <- data.frame(gene_id = rownames(dds)) %>%
  left_join(mapping, by = "gene_id") %>%
  filter(!is.na(Entrez)) %>%
  mutate(Entrez = as.character(Entrez)) %>%
  pull(Entrez) %>%
  unique()

message("Background size: ", length(background_genes), " genes")

# ============================================================
# 7) Load query gene list (common DEGs WW)
# ============================================================
gene_list <- read_excel("DEG_results/effetto_geno/WW/common_deg_genes_WW.xlsx") %>%
  as.data.frame()

deg_mapped <- gene_list %>% left_join(mapping, by = "gene_id")

query <- deg_mapped %>%
  filter(!is.na(Entrez)) %>%
  pull(Entrez) %>%
  unique()

if (length(query) < 5) {
  stop("Too few genes mapped to Entrez (", length(query), "). Minimum 5 required.")
}

message("Query size: ", length(query), " genes")

# ============================================================
# 8) g:Profiler enrichment
# ============================================================
outdir <- "GO_enrichment_results_WW"
dir.create(outdir, showWarnings = FALSE)

res_gp <- gost(
  query             = query,
  organism          = "slycopersicum",
  user_threshold    = 0.05,
  correction_method = "fdr",
  domain_scope      = "custom_annotated",
  custom_bg         = background_genes
)

if (is.null(res_gp$result)) stop("gost() returned no significant results.")

res_df <- as.data.frame(res_gp$result)
write_xlsx(res_df, file.path(outdir, "GO_enrichment_full_common_genes.xlsx"))

# ============================================================
# 9) Per-ontology output + rrvgo clustering + barplots
# ============================================================
ontologies <- c("GO:BP", "GO:MF", "GO:CC", "KEGG")

for (ont in ontologies) {

  res_sub  <- res_df %>% filter(source == ont)
  if (nrow(res_sub) == 0) next

  safe_ont <- gsub("[:]", "_", ont)

  write_xlsx(res_sub,
             file.path(outdir, paste0("GO_enrichment_", safe_ont, "_common_genes.xlsx")))

  # rrvgo redundancy reduction — GO ontologies only, minimum 5 terms
  if (grepl("GO:", ont) && nrow(res_sub) > 5) {

    message("Running rrvgo for: ", ont)

    simMatrix <- calculateSimMatrix(
      res_sub$term_id,
      orgdb  = org.Slycopersicum.eg.db,
      ont    = sub("GO:", "", ont),
      method = "Rel"
    )

    scores  <- setNames(-log10(res_sub$p_value), res_sub$term_id)
    reduced <- reduceSimMatrix(simMatrix, scores,
                               threshold = 0.7,
                               orgdb     = org.Slycopersicum.eg.db)
    reduced$parents <- NULL

    write_xlsx(reduced,
               file.path(outdir, paste0("GO_enrichment_reduced_", safe_ont, "_common_genes.xlsx")))

    pdf(file.path(outdir, paste0("GO_treemap_", safe_ont, "_common_genes.pdf")),
        width = 8, height = 6)
    treemapPlot(reduced)
    dev.off()
  }

  # Barplot — top 20 terms by FDR
  top_terms <- res_sub %>% arrange(p_value) %>% head(20)

  p <- ggplot(top_terms,
              aes(x = reorder(term_name, -log10(p_value)),
                  y = -log10(p_value))) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    theme_bw() +
    labs(title = paste("Top terms —", ont, "— common DEGs WW"),
         x = "Term", y = "-log10(FDR)")

  ggsave(file.path(outdir, paste0("GO_top_terms_", safe_ont, "_common_genes.pdf")),
         p, width = 8, height = 6)
}
