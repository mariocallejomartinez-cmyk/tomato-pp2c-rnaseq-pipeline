# ============================================================
# Build org.Slycopersicum.eg.db annotation package from NCBI
# Run once before the enrichment analysis (04_gprofiler_rrvgo.R)
# Requires internet connection — downloads gene data from NCBI
# Tax ID 4081 = Solanum lycopersicum
#
# Known issue: makeOrgPackageFromNCBI() may fail on macOS due
# to URL handling differences. Run this script on Windows if
# the macOS build fails.
# ============================================================

library(AnnotationForge)

makeOrgPackageFromNCBI(
  version    = "1.0",
  maintainer = "Your Name <your.email@institution.it>",
  author     = "Your Name <your.email@institution.it>",
  outputDir  = ".",
  tax_id     = "4081",
  genus      = "Solanum",
  species    = "lycopersicum"
)

# Install the generated package locally
install.packages("./org.Slycopersicum.eg.db", repos = NULL, type = "source")

# The package can then be loaded in other scripts with:
# library(org.Slycopersicum.eg.db)
