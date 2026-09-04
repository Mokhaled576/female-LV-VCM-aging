# ============================================================
# run_all.R
# Master driver for the computational analysis
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n[1/4] Defining female LV cohort...\n")
source("scripts/01_metadata_cohort.R")

cat("\n[2/4] Preparing/subsetting expression data...\n")
source("scripts/02_load_expression_data.R")

cat("\n[3/4] QC, normalization, PCA, UMAP and Harmony...\n")
source("scripts/03_VCM_QC_preprocessing.R")

cat("\n[4/4] Donor pseudobulk DE, sensitivity and enrichment...\n")
source("scripts/04_VCM_pseudobulk_age_DE.R")

cat("\nWorkflow complete.\n")
