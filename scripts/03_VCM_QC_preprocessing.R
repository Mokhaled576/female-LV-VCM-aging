# ============================================================
# 03_VCM_QC_preprocessing.R
# GitHub-ready version: project-relative paths
# ============================================================

# ============================================================
# Project: Healthy female human heart aging analysis
# Script 03: Ventricular cardiomyocyte QC and preprocessing
#
# Dataset:
# Read et al. 2024 human heart snRNA-seq meta-analysis
#
# Primary population:
# Female donors
# Left ventricle
# Ventricular cardiomyocytes
#
# Primary biological comparison:
# Younger_20_45 vs Older_55plus
#
# IMPORTANT:
# This script follows the course workflow:
#
# 1. Build single-cell object
# 2. Per-nucleus QC
# 3. QC visualization
# 4. Define justified QC thresholds
# 5. Filter nuclei if needed
# 6. Normalize
# 7. Identify highly variable genes
# 8. PCA
# 9. UMAP
# 10. Assess donor/study effects
# 11. Decide whether integration is needed
#
# Differential expression is NOT performed here.
# Final DE will use raw counts aggregated to donor-level
# pseudobulk samples.
# ============================================================


# ------------------------------------------------------------
# STEP 1: Load required packages
# ------------------------------------------------------------

library(SingleCellExperiment)
library(scuttle)
library(Matrix)
library(edgeR)

# Load tidyverse last because some Bioconductor packages
# can mask commonly used functions.
library(tidyverse)


# ------------------------------------------------------------
# STEP 2: Set a random seed
# ------------------------------------------------------------
# PCA itself is deterministic, but later UMAP and some other
# dimensionality-reduction procedures can be stochastic.
#
# We record the seed for reproducibility.

set.seed(12345)


# ------------------------------------------------------------
# STEP 3: Confirm project directory
# ------------------------------------------------------------

getwd()


# ------------------------------------------------------------
# STEP 4: Define input files
# ------------------------------------------------------------
# These files were created or downloaded during Script 02.

count_file <- "data/female_LV_VCM_counts.mtx"

cell_meta_file <- "results/tables/female_LV_VCM_cell_metadata.csv"

gene_meta_file <- file.path("data_raw", "multi_study_snRNA_Seq_gene_data.csv")


# ------------------------------------------------------------
# STEP 5: Confirm all required files exist
# ------------------------------------------------------------

file.exists(count_file)

file.exists(cell_meta_file)

file.exists(gene_meta_file)
# ------------------------------------------------------------
# STEP 6: Load target-cell metadata
# ------------------------------------------------------------
# One row = one LV ventricular cardiomyocyte nucleus.
#
# Expected:
# 72,104 rows

target_cell_meta <- readr::read_csv(
  cell_meta_file,
  show_col_types = FALSE
)

dim(target_cell_meta)

names(target_cell_meta)
# ------------------------------------------------------------
# STEP 7: Load gene metadata
# ------------------------------------------------------------

gene_meta <- readr::read_csv(
  gene_meta_file,
  show_col_types = FALSE
)

dim(gene_meta)

names(gene_meta)

# ------------------------------------------------------------
# STEP 8: Create unique gene identifiers
# ------------------------------------------------------------
# gene_short_name contains one duplicated symbol:
# TMSB15B
#
# We preserve both rows and create unique computational
# row names using make.unique().

gene_meta <- gene_meta %>%
  dplyr::mutate(
    gene_id_unique = make.unique(gene_short_name)
  )


# Confirm uniqueness

sum(
  duplicated(gene_meta$gene_id_unique)
)

dplyr::n_distinct(
  gene_meta$gene_id_unique
)
# ------------------------------------------------------------
# STEP 9: Load raw VCM count matrix
# ------------------------------------------------------------
# Matrix dimensions:
#
# 17,726 genes x 72,104 nuclei
#
# IMPORTANT:
# Keep this matrix sparse.
# Never convert it with as.matrix().

vcm_counts <- Matrix::readMM(
  count_file
)


# ------------------------------------------------------------
# STEP 10: Confirm matrix dimensions
# ------------------------------------------------------------

dim(vcm_counts)

class(vcm_counts)

# ------------------------------------------------------------
# STEP 11: Attach unique gene identifiers
# ------------------------------------------------------------

rownames(vcm_counts) <- gene_meta$gene_id_unique


# Confirm gene alignment

stopifnot(
  nrow(vcm_counts) == nrow(gene_meta)
)

stopifnot(
  identical(
    rownames(vcm_counts),
    gene_meta$gene_id_unique
  )
)
# ------------------------------------------------------------
# STEP 12: Confirm cell metadata alignment
# ------------------------------------------------------------

stopifnot(
  ncol(vcm_counts) == nrow(target_cell_meta)
)

cat(
  "SUCCESS: raw counts align with target-cell metadata.\n"
)

# ------------------------------------------------------------
# STEP 13: Create SingleCellExperiment object
# ------------------------------------------------------------
#
# assays:
# raw counts
#
# colData:
# nucleus-level metadata
#
# rowData:
# original and unique gene identifiers

sce_vcm <- SingleCellExperiment::SingleCellExperiment(
  assays = list(
    counts = vcm_counts
  ),
  colData = S4Vectors::DataFrame(
    target_cell_meta
  ),
  rowData = S4Vectors::DataFrame(
    gene_meta
  )
)
# ------------------------------------------------------------
# STEP 14: Inspect the SingleCellExperiment
# ------------------------------------------------------------

sce_vcm

dim(sce_vcm)

assayNames(sce_vcm)

colnames(colData(sce_vcm))

colnames(rowData(sce_vcm))

# ------------------------------------------------------------
# STEP 15: Identify mitochondrial genes
# ------------------------------------------------------------
# Human mitochondrial genes are typically named with
# the prefix "MT-".
#
# We use the original biological gene symbol column,
# not the make.unique() computational ID.

is_mito <- stringr::str_detect(
  rowData(sce_vcm)$gene_short_name,
  "^MT-"
)


# Number of mitochondrial genes represented

sum(is_mito)


# Inspect mitochondrial gene names

rowData(sce_vcm)$gene_short_name[
  is_mito
]
# ------------------------------------------------------------
# STEP 16: Calculate per-nucleus QC metrics
# ------------------------------------------------------------
#
# This calculates:
#
# sum:
# total UMI counts per nucleus
#
# detected:
# number of genes detected per nucleus
#
# subsets_Mito_percent:
# percentage of counts from mitochondrial genes

qc_metrics <- scuttle::perCellQCMetrics(
  sce_vcm,
  subsets = list(
    Mito = is_mito
  )
)
dim(qc_metrics)

# ------------------------------------------------------------
# STEP 17: Attach QC metrics to colData
# ------------------------------------------------------------

colData(sce_vcm)$total_counts <- qc_metrics$sum

colData(sce_vcm)$detected_genes <- qc_metrics$detected

colData(sce_vcm)$mito_percent <- qc_metrics$subsets_Mito_percent
# ------------------------------------------------------------
# STEP 18: Summarize per-nucleus QC metrics
# ------------------------------------------------------------

summary(
  colData(sce_vcm)$total_counts
)

summary(
  colData(sce_vcm)$detected_genes
)

summary(
  colData(sce_vcm)$mito_percent
)
# ------------------------------------------------------------
# STEP 19: Summarize QC metrics by donor
# ------------------------------------------------------------

qc_by_donor <- as.data.frame(
  colData(sce_vcm)
) %>%
  dplyr::group_by(
    Donor,
    Age,
    DataSource
  ) %>%
  dplyr::summarise(
    n_nuclei = dplyr::n(),
    
    
    # ------------------------------------------------------------
    # STEP 20: Add age-group labels to the nucleus metadata
    # ------------------------------------------------------------
    # The deposited cell metadata contains continuous age but
    # does not contain our custom age_group variable.
    #
    # We recreate exactly the same primary groups:
    #
    # Younger: 20-45 years
    # Older:   >=55 years
    #
    # There should be no 46-54 year nuclei because those donors
    # were excluded when the primary cohort was constructed.
    
    colData(sce_vcm)$age_group <- dplyr::case_when(
      
      colData(sce_vcm)$Age >= 20 &
        colData(sce_vcm)$Age <= 45 ~ "Younger_20_45",
      
      colData(sce_vcm)$Age >= 55 ~ "Older_55plus",
      
      TRUE ~ NA_character_)
    
    
    # Confirm that every nucleus has an age-group assignment.
    
    table(
      colData(sce_vcm)$age_group,
      useNA = "ifany"
      )
    
    median_counts = median(total_counts),
    median_detected_genes = median(detected_genes),
    median_mito_percent = median(mito_percent),
    
    min_counts = min(total_counts),
    max_counts = max(total_counts),
    
    min_detected_genes = min(detected_genes),
    max_detected_genes = max(detected_genes),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    median_counts
  )

print(
  qc_by_donor,
  n = Inf,
  width = Inf
)
# ------------------------------------------------------------
# STEP 21: Create QC plotting data frame
# ------------------------------------------------------------

qc_plot_data <- as.data.frame(
  colData(sce_vcm)
)
dim(qc_plot_data)
# ------------------------------------------------------------
# STEP 22: Plot total counts per nucleus
# ------------------------------------------------------------

p_counts <- ggplot(
  qc_plot_data,
  aes(
    x = total_counts
  )
) +
  geom_histogram(
    bins = 100
  ) +
  scale_x_log10() +
  labs(
    title = "VCM nuclei: total UMI counts",
    x = "Total counts per nucleus (log10 scale)",
    y = "Number of nuclei"
  ) +
  theme_classic()

p_counts

# ------------------------------------------------------------
# STEP 23: Plot detected genes per nucleus
# ------------------------------------------------------------

p_genes <- ggplot(
  qc_plot_data,
  aes(
    x = detected_genes
  )
) +
  geom_histogram(
    bins = 100
  ) +
  scale_x_log10() +
  labs(
    title = "VCM nuclei: detected genes",
    x = "Detected genes per nucleus (log10 scale)",
    y = "Number of nuclei"
  ) +
  theme_classic()

p_genes

# ------------------------------------------------------------
# STEP 24: Plot mitochondrial percentage
# ------------------------------------------------------------

p_mito <- ggplot(
  qc_plot_data,
  aes(
    x = mito_percent
  )
) +
  geom_histogram(
    bins = 100
  ) +
  labs(
    title = "VCM nuclei: mitochondrial percentage",
    x = "Mitochondrial counts (%)",
    y = "Number of nuclei"
  ) +
  theme_classic()

p_mito

# ------------------------------------------------------------
# STEP 25: Compare detected genes across donors
# ------------------------------------------------------------
# Donors are ordered by median detected genes.
#
# This lets us see whether a global threshold would
# disproportionately remove nuclei from particular donors.

donor_order_genes <- qc_plot_data %>%
  dplyr::group_by(
    Donor
  ) %>%
  dplyr::summarise(
    median_detected = median(detected_genes),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    median_detected
  ) %>%
  dplyr::pull(
    Donor
  )

qc_plot_data$Donor_genes_ordered <- factor(
  qc_plot_data$Donor,
  levels = donor_order_genes
)


p_genes_donor <- ggplot(
  qc_plot_data,
  aes(
    x = Donor_genes_ordered,
    y = detected_genes
  )
) +
  geom_boxplot(
    outlier.size = 0.2
  ) +
  scale_y_log10() +
  labs(
    title = "Detected genes per VCM nucleus by donor",
    x = "Donor",
    y = "Detected genes (log10 scale)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1
    )
  )

p_genes_donor
# ------------------------------------------------------------
# STEP 26: Compare total counts across donors
# ------------------------------------------------------------

p_counts_donor <- ggplot(
  qc_plot_data,
  aes(
    x = Donor_genes_ordered,
    y = total_counts
  )
) +
  geom_boxplot(
    outlier.size = 0.2
  ) +
  scale_y_log10() +
  labs(
    title = "Total counts per VCM nucleus by donor",
    
    
    x = "Donor",
    y = "Total counts (log10 scale)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1
    )
  )

p_counts_donor
# ------------------------------------------------------------
# STEP 27: Compare detected genes across source studies
# ------------------------------------------------------------

p_genes_study <- ggplot(
  qc_plot_data,
  aes(
    x = DataSource,
    y = detected_genes
  )
) +
  geom_boxplot(
    outlier.size = 0.2
  ) +
  scale_y_log10() +
  labs(
    title = "Detected genes per VCM nucleus by study",
    x = "Original study",
    y = "Detected genes (log10 scale)"
  ) +
  theme_classic()

p_genes_study
# ------------------------------------------------------------
# STEP 28: Compare total counts across source studies
# ------------------------------------------------------------

p_counts_study <- ggplot(
  qc_plot_data,
  aes(
    x = DataSource,
    y = total_counts
  )
) +
  geom_boxplot(
    outlier.size = 0.2
  ) +
  scale_y_log10() +
  labs(
    title = "Total counts per VCM nucleus by study",
    x = "Original study",
    y = "Total counts (log10 scale)"
  ) +
  theme_classic()

p_counts_study

# ------------------------------------------------------------
# STEP 29: Identify unusually low-count nuclei within donor
# ------------------------------------------------------------
# nmads = 3:
# nucleus must lie more than 3 median absolute deviations
# below its donor-specific distribution.
#
# log = TRUE:
# library size is evaluated on a log scale because its
# distribution is strongly right-skewed.
#
# batch = Donor:
# thresholds are estimated separately for each donor.
#
# IMPORTANT:
# This only FLAGS potential low-quality nuclei.
# Nothing is removed yet.

low_count_flag <- scuttle::isOutlier(
  qc_plot_data$total_counts,
  type = "lower",
  log = TRUE,
  nmads = 3,
  batch = qc_plot_data$Donor
)
# ------------------------------------------------------------
# STEP 30: Identify unusually low-feature nuclei within donor
# ------------------------------------------------------------

low_gene_flag <- scuttle::isOutlier(
  qc_plot_data$detected_genes,
  type = "lower",
  log = TRUE,
  nmads = 3,
  batch = qc_plot_data$Donor
)
# ------------------------------------------------------------
# STEP 31: Identify mitochondrial-percentage outliers
# ------------------------------------------------------------
# Again, this is donor-aware.
#
# We are not imposing an arbitrary 5%, 10%, or 20% cutoff.

high_mito_flag <- scuttle::isOutlier(
  qc_plot_data$mito_percent,
  type = "higher",
  nmads = 3,
  batch = qc_plot_data$Donor
)
# ------------------------------------------------------------
# STEP 32: Count QC flags
# ------------------------------------------------------------

sum(low_count_flag)

sum(low_gene_flag)

sum(high_mito_flag)
# ------------------------------------------------------------
# STEP 33: Create preliminary combined QC flag
# ------------------------------------------------------------
#
# A nucleus is provisionally flagged if it has:
#
# - unusually low counts OR
# - unusually few genes OR
# - unusually high mitochondrial percentage
#
# Still no filtering yet.

qc_fail_preliminary <- (
  low_count_flag |
    low_gene_flag |
    high_mito_flag
)

sum(qc_fail_preliminary)

mean(qc_fail_preliminary) * 100

# ------------------------------------------------------------
# STEP 34: Examine QC failures by donor
# ------------------------------------------------------------

qc_flag_table <- qc_plot_data %>%
  dplyr::mutate(
    low_counts = low_count_flag,
    low_genes = low_gene_flag,
    high_mito = high_mito_flag,
    qc_fail = qc_fail_preliminary
  ) %>%
  dplyr::group_by(
    Donor,
    Age,
    DataSource,
    age_group
  ) %>%
  dplyr::summarise(
    
    n_before = dplyr::n(),
    
    n_low_counts = sum(low_counts),
    
    n_low_genes = sum(low_genes),
    
    n_high_mito = sum(high_mito),
    
    n_qc_fail = sum(qc_fail),
    
    n_pass = sum(!qc_fail),
    
    percent_fail = 100 * mean(qc_fail),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(percent_fail)
  )


print(
  qc_flag_table,
  n = Inf,
  width = Inf
)
# ------------------------------------------------------------
# STEP 35: Examine mitochondrial percentage by donor
# ------------------------------------------------------------
# The preliminary MAD-based mitochondrial rule flagged many
# nuclei even though mitochondrial percentages are extremely
# low overall.
#
# We therefore inspect the actual mitochondrial distributions
# before using mitochondrial percentage as a filtering rule.

mito_by_donor <- qc_plot_data %>%
  dplyr::group_by(
    Donor,
    Age,
    DataSource,
    age_group
  ) %>%
  dplyr::summarise(
    
    n_nuclei = dplyr::n(),
    
    median_mito = median(mito_percent),
    
    p75_mito = quantile(
      mito_percent,
      probs = 0.75
    ),
    
    p90_mito = quantile(
      mito_percent,
      probs = 0.90
    ),
    
    p95_mito = quantile(
      mito_percent,
      probs = 0.95
    ),
    
    p99_mito = quantile(
      mito_percent,
      probs = 0.99
    ),
    
    max_mito = max(mito_percent),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(p95_mito)
  )


print(
  mito_by_donor,
  n = Inf,
  width = Inf
)
# ------------------------------------------------------------
# STEP 36: Count nuclei above several mitochondrial percentages
# ------------------------------------------------------------
# These values help us understand the scale of the
# mitochondrial distribution.
#
# They are NOT filtering thresholds.

mito_diagnostic <- tibble::tibble(
  
  cutoff_percent = c(
    0.5,
    1,
    2,
    5,
    10
  ),
  
  n_nuclei = c(
    sum(qc_plot_data$mito_percent > 0.5),
    sum(qc_plot_data$mito_percent > 1),
    sum(qc_plot_data$mito_percent > 2),
    sum(qc_plot_data$mito_percent > 5),
    sum(qc_plot_data$mito_percent > 10)
  )
  
) %>%
  dplyr::mutate(
    percent_nuclei =
      100 * n_nuclei / nrow(qc_plot_data)
  )


mito_diagnostic

# ------------------------------------------------------------
# STEP 37: Examine counts versus detected genes
# ------------------------------------------------------------

p_counts_vs_genes <- ggplot(
  qc_plot_data,
  aes(
    x = total_counts,
    y = detected_genes
  )
) +
  geom_point(
    alpha = 0.15,
    size = 0.3
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "VCM nucleus complexity",
    x = "Total counts per nucleus (log10)",
    y = "Detected genes per nucleus (log10)"
  ) +
  theme_classic()

p_counts_vs_genes

# ------------------------------------------------------------
# STEP 38: Evaluate low-count and low-gene QC independently
# ------------------------------------------------------------

low_complexity_flag <- (
  low_count_flag |
    low_gene_flag
)


sum(low_complexity_flag)

100 * mean(low_complexity_flag)

# ------------------------------------------------------------
# STEP 39: Low-complexity flags by donor
# ------------------------------------------------------------

low_complexity_table <- qc_plot_data %>%
  dplyr::mutate(
    low_count = low_count_flag,
    low_gene = low_gene_flag,
    low_complexity = low_complexity_flag
  ) %>%
  dplyr::group_by(
    Donor,
    Age,
    DataSource,
    age_group
  ) %>%
  dplyr::summarise(
    
    n_before = dplyr::n(),
    
    n_low_count = sum(low_count),
    
    n_low_gene = sum(low_gene),
    
    n_low_complexity = sum(low_complexity),
    
    n_remaining =
      n_before - n_low_complexity,
    
    percent_low_complexity =
      100 * n_low_complexity / n_before,
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(percent_low_complexity)
  )


print(
  low_complexity_table,
  n = Inf,
  width = Inf
)
qc_pass <- !low_complexity_flag
# ------------------------------------------------------------
# STEP 40: Define the final secondary QC rule
# ------------------------------------------------------------
#
# FINAL DECISION:
#
# A nucleus fails secondary QC if it is an unusually
# low-complexity nucleus within its own donor.
#
# Specifically, it is removed if it is:
#
# 1. a donor-specific low-count outlier, OR
# 2. a donor-specific low-detected-gene outlier
#
# Both were defined using:
#
# - 3 median absolute deviations
# - lower tail
# - log-transformed metric
# - thresholds calculated separately within donor
#
# Mitochondrial percentage is NOT used as an exclusion
# criterion because mitochondrial fractions are extremely
# low in this snRNA-seq dataset and donor-specific MAD
# thresholds falsely flagged large fractions of otherwise
# plausible nuclei.

qc_pass <- !low_complexity_flag
# ------------------------------------------------------------

# ------------------------------------------------------------
# STEP 42: Report nuclei before and after QC
# ------------------------------------------------------------

n_before_qc <- ncol(sce_vcm)

n_pass_qc <- sum(
  colData(sce_vcm)$qc_pass
)

n_fail_qc <- sum(
  !colData(sce_vcm)$qc_pass
)

percent_removed_qc <- 100 * n_fail_qc / n_before_qc


cat(
  "Nuclei before QC:", n_before_qc, "\n"
)

cat(
  "Nuclei retained after QC:", n_pass_qc, "\n"
)

cat(
  "Nuclei removed:", n_fail_qc, "\n"
)

cat(
  "Percent removed:", round(percent_removed_qc, 2), "%\n"
)

# ------------------------------------------------------------
# STEP 43: QC retention by age group
# ------------------------------------------------------------

qc_age_table <- as.data.frame(
  colData(sce_vcm)
) %>%
  dplyr::group_by(
    age_group
  ) %>%
  dplyr::summarise(
    
    n_before = dplyr::n(),
    
    n_pass = sum(qc_pass),
    
    n_removed = sum(!qc_pass),
    
    percent_removed =
      100 * mean(!qc_pass),
    
    .groups = "drop"
  )


qc_age_table
# ------------------------------------------------------------
# STEP 44: Create QC-filtered SCE
# ------------------------------------------------------------

sce_vcm_qc <- sce_vcm[
  ,
  colData(sce_vcm)$qc_pass
]
dim(sce_vcm_qc)
# ------------------------------------------------------------
# STEP 45: Confirm all donors remain after QC
# ------------------------------------------------------------

length(
  unique(
    colData(sce_vcm_qc)$Donor
  )
)

table(
  colData(sce_vcm_qc)$Donor
)

# ------------------------------------------------------------
# STEP 46: Save QC summary tables
# ------------------------------------------------------------

readr::write_csv(
  low_complexity_table,
  "results/tables/VCM_QC_by_donor.csv"
)

readr::write_csv(
  qc_age_table,
  "results/tables/VCM_QC_by_age_group.csv"
)
# ------------------------------------------------------------
# STEP 47: Save metadata for QC-passing nuclei
# ------------------------------------------------------------

qc_pass_metadata <- as.data.frame(
  colData(sce_vcm_qc)
)

readr::write_csv(
  qc_pass_metadata,
  "results/tables/female_LV_VCM_QC_pass_metadata.csv"
)
# ------------------------------------------------------------
# STEP 48: Final QC sanity check
# ------------------------------------------------------------

cat(
  "Genes:", nrow(sce_vcm_qc), "\n"
)

cat(
  "QC-passing nuclei:", ncol(sce_vcm_qc), "\n"
)

cat(
  "Donors:",
  length(unique(colData(sce_vcm_qc)$Donor)),
  "\n"
)

cat(
  "Younger nuclei:",
  sum(colData(sce_vcm_qc)$age_group == "Younger_20_45"),
  "\n"
)

cat(
  "Older nuclei:",
  sum(colData(sce_vcm_qc)$age_group == "Older_55plus"),
  "\n"
)
# ------------------------------------------------------------
# STEP 49: Remove genes with zero counts after QC
# ------------------------------------------------------------
# Some genes could theoretically become completely absent
# after removing low-quality nuclei.
#
# We remove only genes with zero total counts across all
# QC-passing nuclei.

gene_totals_qc <- Matrix::rowSums(
  counts(sce_vcm_qc)
)

sum(gene_totals_qc == 0)
# Keep genes detected at least once.

sce_vcm_qc <- sce_vcm_qc[
  gene_totals_qc > 0,
]

dim(sce_vcm_qc)
# ------------------------------------------------------------
# STEP 50: Compute library-size factors
# ------------------------------------------------------------
#
# For exploratory normalization, start with library-size
# normalization.
#
# The size factor represents the relative sequencing depth
# of each nucleus.

sce_vcm_qc <- scuttle::computeLibraryFactors(
  sce_vcm_qc
)
summary(
  sizeFactors(sce_vcm_qc)
)
# ------------------------------------------------------------
# STEP 51: Log-normalize expression
# ------------------------------------------------------------
#
# logNormCounts() produces normalized log-expression values.
#
# Raw counts remain preserved in the "counts" assay.
# The normalized values are stored separately in "logcounts".
#
# This is important:
# edgeR pseudobulk DE later will use RAW counts,
# not these normalized values.

sce_vcm_qc <- scuttle::logNormCounts(
  sce_vcm_qc
)
assayNames(sce_vcm_qc)
# ------------------------------------------------------------
# STEP 52: Load scran
# ------------------------------------------------------------

library(scran)
# ------------------------------------------------------------
# STEP 53: Model gene variance
# ------------------------------------------------------------
#
# modelGeneVar() separates biological variation from
# technical variation in normalized log-expression.
#
# We will use this to choose highly variable genes for PCA.

gene_variance <- scran::modelGeneVar(
  sce_vcm_qc
)
gene_variance
summary(
  gene_variance$bio
)
# ------------------------------------------------------------
# STEP 54: Select top highly variable genes
# ------------------------------------------------------------
#
# Use the 2,000 genes with the strongest modeled
# biological variation.
#
# This keeps PCA focused on informative genes and limits
# memory usage.

n_hvg <- min(
  2000,
  nrow(gene_variance)
)

hvg_order <- order(
  gene_variance$bio,
  decreasing = TRUE
)

hvg_genes <- rownames(
  gene_variance
)[
  hvg_order[
    seq_len(n_hvg)
  ]
]


length(hvg_genes)

head(
  hvg_genes,
  20
)
# ------------------------------------------------------------
# STEP 55: Check biologically relevant genes among HVGs
# ------------------------------------------------------------
#
# PURPOSE:
#
# Our primary biological interests include:
#
# 1. Cytochrome P450 enzymes
# 2. EET / epoxyeicosatrienoic acid pathway
# 3. HETE / hydroxyeicosatetraenoic acid pathway
# 4. Soluble epoxide hydrolase pathway
# 5. Arachidonic acid / eicosanoid metabolism
# 6. Estrogen receptors, synthesis, and metabolism
# 7. PPAR signaling
# 8. AHR signaling
# 9. Xenobiotic/nuclear receptor signaling
# 10. Oxidative stress and related metabolic regulators
#
# IMPORTANT:
#
# This is ONLY a diagnostic check.
#
# We do NOT force any of these genes into the HVG list.
# HVGs remain selected completely independently from the
# biological hypothesis.
# ------------------------------------------------------------


# ------------------------------------------------------------
# STEP 55A: Capture ALL CYP genes present in the dataset
# ------------------------------------------------------------
#
# Rather than manually choosing only a few CYP genes,
# identify every gene beginning with "CYP" that occurs
# in the Read et al. expression matrix.

all_cyp_genes_in_dataset <- rowData(sce_vcm_qc)$gene_short_name[
  stringr::str_detect(
    rowData(sce_vcm_qc)$gene_short_name,
    "^CYP"
  )
]

all_cyp_genes_in_dataset <- sort(
  unique(all_cyp_genes_in_dataset)
)


# See which CYP genes are available

all_cyp_genes_in_dataset
# ------------------------------------------------------------
# STEP 55B: EET / epoxygenase / epoxide-hydrolase pathway
# ------------------------------------------------------------
#
# CYP2J2 and several CYP2C enzymes can generate EETs from
# arachidonic acid.
#
# EPHX2 encodes soluble epoxide hydrolase (sEH), which
# converts EETs to less-active DHETs.
#
# EPHX1 encodes microsomal epoxide hydrolase.

eet_epoxide_genes <- c(
  
  "CYP2J2",
  
  "CYP2C8",
  "CYP2C9",
  "CYP2C18",
  "CYP2C19",
  
  "EPHX2",       # soluble epoxide hydrolase / sEH
  "EPHX1",       # microsomal epoxide hydrolase
  
  "POR",         # NADPH-cytochrome P450 reductase
  "CYB5A",       # cytochrome b5
  "CYB5B"
)
# ------------------------------------------------------------
# STEP 55C: HETE / CYP omega-hydroxylase pathway
# ------------------------------------------------------------
#
# Includes CYP enzymes associated with formation and/or
# metabolism of HETEs such as 20-HETE.

hete_cyp_genes <- c(
  
  "CYP4A11",
  "CYP4A22",
  
  "CYP4F2",
  "CYP4F3",
  "CYP4F8",
  "CYP4F11",
  "CYP4F12",
  "CYP4F22",
  
  "CYP2U1",
  "CYP2J2",
  
  "POR",
  "CYB5A",
  
  # Proposed/recognized signaling receptor associated
  # with 20-HETE signaling
  "GPR75"
)
# ------------------------------------------------------------
# STEP 55D: Arachidonic acid release and metabolism
# ------------------------------------------------------------

arachidonic_acid_genes <- c(
  
  # ----------------------------------------------------------
  # Phospholipases releasing arachidonic acid
  # ----------------------------------------------------------
  
  "PLA2G4A",
  "PLA2G4B",
  "PLA2G4C",
  "PLA2G4D",
  "PLA2G4E",
  "PLA2G4F",
  
  "PLA2G6",
  "PLA2G2A",
  "PLA2G2D",
  "PLA2G2E",
  "PLA2G2F",
  
  
  # ----------------------------------------------------------
  # Cyclooxygenase pathway
  # ----------------------------------------------------------
  
  "PTGS1",
  "PTGS2",
  
  
  # ----------------------------------------------------------
  # Prostaglandin synthases
  # ----------------------------------------------------------
  
  "PTGES",
  "PTGES2",
  "PTGES3",
  
  "PTGDS",
  "HPGDS",
  
  "PTGIS",
  
  "TBXAS1",
  
  
  # ----------------------------------------------------------
  # Prostaglandin degradation
  # ----------------------------------------------------------
  
  "HPGD",
  "PTGR1",
  "PTGR2"
)
# ------------------------------------------------------------
# STEP 55E: Lipoxygenase / HETE / leukotriene pathway
# ------------------------------------------------------------

lox_hete_genes <- c(
  
  "ALOX5",
  "ALOX5AP",
  
  "ALOX12",
  "ALOX12B",
  
  "ALOX15",
  "ALOX15B",
  
  "LTA4H",
  
  "LTC4S",
  
  "MGST1",
  "MGST2",
  "MGST3",
  
  "GGT1",
  "GGT5"
)
# ------------------------------------------------------------
# STEP 55F: Eicosanoid receptors
# ------------------------------------------------------------

eicosanoid_receptor_genes <- c(
  
  # PGE receptors
  "PTGER1",
  "PTGER2",
  "PTGER3",
  "PTGER4",
  
  # PGD receptor
  "PTGDR",
  "PTGDR2",
  
  # Prostacyclin receptor
  "PTGIR",
  
  # Thromboxane receptor
  "TBXA2R",
  
  # Leukotriene B4 receptors
  "LTB4R",
  "LTB4R2",
  
  # Cysteinyl leukotriene receptors
  "CYSLTR1",
  "CYSLTR2",
  
  # Oxoeicosanoid receptor
  "OXER1",
  
  # 20-HETE-associated receptor
  "GPR75"
)
# ------------------------------------------------------------
# STEP 55G: Estrogen receptor/signaling genes
# ------------------------------------------------------------

estrogen_signaling_genes <- c(
  
  # Classical nuclear estrogen receptors
  "ESR1",
  "ESR2",
  
  # Membrane estrogen receptor
  "GPER1",
  
  # Estrogen-related receptors
  "ESRRA",
  "ESRRB",
  "ESRRG",
  
  # Nuclear receptor coactivators/corepressors
  "NCOA1",
  "NCOA2",
  "NCOA3",
  
  "NCOR1",
  "NCOR2"
)
# ------------------------------------------------------------
# STEP 55H: Estrogen synthesis and metabolism
# ------------------------------------------------------------

estrogen_metabolism_genes <- c(
  
  # Aromatase
  "CYP19A1",
  
  # Steroidogenic CYP enzymes
  "CYP11A1",
  "CYP17A1",
  "CYP21A2",
  
  # Hydroxysteroid dehydrogenases
  "HSD17B1",
  "HSD17B2",
  "HSD17B4",
  "HSD17B7",
  "HSD17B10",
  "HSD17B11",
  "HSD17B12",
  "HSD17B13",
  "HSD17B14",
  
  "HSD3B1",
  "HSD3B2",
  
  # Steroid sulfatase
  "STS",
  
  # Estrogen sulfotransferase
  "SULT1E1",
  
  # Catechol-estrogen metabolism
  "COMT",
  
  # Relevant glucuronidation enzymes
  "UGT1A1",
  "UGT1A3",
  "UGT1A4",
  "UGT1A6",
  "UGT1A7",
  "UGT1A8",
  "UGT1A9",
  "UGT1A10",
  
  "UGT2B4",
  "UGT2B7",
  "UGT2B10",
  "UGT2B15",
  "UGT2B17"
)
# ------------------------------------------------------------
# STEP 55I: PPAR / lipid-metabolism signaling
# ------------------------------------------------------------

ppar_genes <- c(
  
  # PPAR receptors
  "PPARA",
  "PPARD",
  "PPARG",
  
  # RXR heterodimer partners
  "RXRA",
  "RXRB",
  "RXRG",
  
  # Major transcriptional coactivators
  "PPARGC1A",
  "PPARGC1B",
  
  # Nuclear receptor coactivators/corepressors
  "NCOA1",
  "NCOA2",
  "NCOA3",
  
  "NCOR1",
  "NCOR2",
  
  # Representative PPAR-regulated cardiac metabolic genes
  "CD36",
  "FABP3",
  "FABP4",
  "CPT1A",
  "CPT1B",
  "CPT2",
  "ACOX1",
  "ACOX2",
  "ACADL",
  "ACADM",
  "ACADVL"
)
# ------------------------------------------------------------
# STEP 55J: AHR signaling
# ------------------------------------------------------------

ahr_genes <- c(
  
  "AHR",
  
  # AHR heterodimer partner
  "ARNT",
  "ARNT2",
  
  # Negative-feedback regulator
  "AHRR",
  
  # Canonical AHR-responsive CYPs
  "CYP1A1",
  "CYP1A2",
  "CYP1B1",
  
  # Other AHR-responsive/regulatory genes
  "TIPARP",
  "NQO1",
  "ALDH3A1"
)
# ------------------------------------------------------------
# STEP 55K: Xenobiotic-sensing nuclear receptors
# ------------------------------------------------------------

xenobiotic_receptor_genes <- c(
  
  # PXR
  "NR1I2",
  
  # CAR
  "NR1I3",
  
  # FXR
  "NR1H4",
  
  # LXR alpha/beta
  "NR1H3",
  "NR1H2",
  
  # HNF4 alpha
  "HNF4A",
  
  # Retinoid X receptors
  "RXRA",
  "RXRB",
  "RXRG"
)
# ------------------------------------------------------------
# STEP 55L: Oxidative stress / redox regulators
# ------------------------------------------------------------

redox_genes <- c(
  
  # NRF2 pathway
  "NFE2L2",
  "KEAP1",
  
  "NQO1",
  
  # Heme oxygenase
  "HMOX1",
  "HMOX2",
  
  # Superoxide dismutases
  "SOD1",
  "SOD2",
  "SOD3",
  
  # Catalase
  "CAT",
  
  # Glutathione peroxidases
  "GPX1",
  "GPX2",
  "GPX3",
  "GPX4",
  
  # Glutathione metabolism
  "GCLC",
  "GCLM",
  "GSR",
  
  # Thioredoxin pathway
  "TXN",
  "TXNRD1",
  "TXNRD2"
)
# ------------------------------------------------------------
# STEP 55M: Combine all genes of interest
# ------------------------------------------------------------

genes_of_interest <- unique(
  c(
    all_cyp_genes_in_dataset,
    
    eet_epoxide_genes,
    
    hete_cyp_genes,
    
    arachidonic_acid_genes,
    
    lox_hete_genes,
    
    eicosanoid_receptor_genes,
    
    estrogen_signaling_genes,
    
    estrogen_metabolism_genes,
    
    ppar_genes,
    
    ahr_genes,
    
    xenobiotic_receptor_genes,
    
    redox_genes
  )
)


length(genes_of_interest)
# ------------------------------------------------------------
# STEP 55N: Check presence in dataset and HVG status
# ------------------------------------------------------------

gene_interest_check <- tibble::tibble(
  
  gene = genes_of_interest,
  
  present_in_dataset =
    genes_of_interest %in%
    rowData(sce_vcm_qc)$gene_short_name,
  
  selected_as_HVG =
    genes_of_interest %in%
    hvg_genes
)


print(
  gene_interest_check,
  n = Inf
)
# ------------------------------------------------------------
# STEP 55O: Determine whether genes are actually detected
#          in the QC-passing VCM nuclei
# ------------------------------------------------------------

gene_total_counts <- Matrix::rowSums(
  counts(sce_vcm_qc)
)


gene_detection_table <- tibble::tibble(
  
  gene = rowData(sce_vcm_qc)$gene_short_name,
  
  total_counts =
    gene_total_counts
  
) %>%
  dplyr::group_by(
    gene
  ) %>%
  dplyr::summarise(
    total_counts = sum(total_counts),
    .groups = "drop"
  )


gene_interest_check <- gene_interest_check %>%
  
  dplyr::left_join(
    gene_detection_table,
    by = "gene"
  ) %>%
  
  dplyr::mutate(
    
    total_counts =
      dplyr::coalesce(
        total_counts,
        0
      ),
    
    detected_in_VCM =
      total_counts > 0
  )


print(
  gene_interest_check,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# Check objects currently available in this R session
# ------------------------------------------------------------

ls()
# ------------------------------------------------------------
# RECOVERY CHECK
# ------------------------------------------------------------

dim(sce_vcm)

assayNames(sce_vcm)

length(unique(colData(sce_vcm)$Donor))

head(
  colnames(colData(sce_vcm))
)
# Check whether the previously calculated QC metrics
# survived inside sce_vcm

c(
  "total_counts",
  "detected_genes",
  "mito_percent",
  "age_group"
) %in% colnames(colData(sce_vcm))
# ------------------------------------------------------------
# RECOVERY STEP 1: Recalculate per-nucleus QC metrics
# ------------------------------------------------------------

# Identify mitochondrial genes from gene symbols
is_mito <- stringr::str_detect(
  rowData(sce_vcm)$gene_short_name,
  "^MT-"
)

# Calculate QC metrics
qc_metrics <- scuttle::perCellQCMetrics(
  sce_vcm,
  subsets = list(Mito = is_mito)
)
# ------------------------------------------------------------
# RECOVERY STEP 3: Recreate donor-specific QC flags
# ------------------------------------------------------------

qc_plot_data <- as.data.frame(
  colData(sce_vcm)
)


low_count_flag <- scuttle::isOutlier(
  qc_plot_data$total_counts,
  type = "lower",
  log = TRUE,
  nmads = 3,
  batch = qc_plot_data$Donor
)


low_gene_flag <- scuttle::isOutlier(
  qc_plot_data$detected_genes,
  type = "lower",
  log = TRUE,
  nmads = 3,
  batch = qc_plot_data$Donor
)


# Mitochondrial percentage is retained as a diagnostic only.
# We are NOT using it to exclude nuclei.

high_mito_flag <- scuttle::isOutlier(
  qc_plot_data$mito_percent,
  type = "higher",
  nmads = 3,
  batch = qc_plot_data$Donor
)


# Final secondary-QC rule:
# remove only low-complexity nuclei

low_complexity_flag <-
  low_count_flag | low_gene_flag

qc_pass <-
  !low_complexity_flag
# Attach the same QC variables used previously
colData(sce_vcm)$total_counts <-
  qc_metrics$sum

colData(sce_vcm)$detected_genes <-
  qc_metrics$detected

colData(sce_vcm)$mito_percent <-
  qc_metrics$subsets_Mito_percent


# ------------------------------------------------------------
# RECOVERY STEP 2: Recreate age groups
# ------------------------------------------------------------

colData(sce_vcm)$age_group <- dplyr::case_when(
  
  colData(sce_vcm)$Age >= 20 &
    colData(sce_vcm)$Age <= 45 ~
    "Younger_20_45",
  
  colData(sce_vcm)$Age >= 55 ~
    "Older_55plus",
  
  TRUE ~ NA_character_
)


# Quick check
table(
  colData(sce_vcm)$age_group,
  useNA = "ifany"
)
# ------------------------------------------------------------
# RECOVERY STEP 4: Verify QC recovery
# ------------------------------------------------------------

cat(
  "Low-count nuclei:",
  sum(low_count_flag),
  "\n"
)

cat(
  "Low-gene nuclei:",
  sum(low_gene_flag),
  "\n"
)

cat(
  "Low-complexity nuclei removed:",
  sum(low_complexity_flag),
  "\n"
)

cat(
  "QC-passing nuclei:",
  sum(qc_pass),
  "\n"
)
# ------------------------------------------------------------
# RECOVERY STEP 5: Recreate final QC-passing SCE
# ------------------------------------------------------------

colData(sce_vcm)$low_count_flag <-
  low_count_flag

colData(sce_vcm)$low_gene_flag <-
  low_gene_flag

colData(sce_vcm)$high_mito_flag <-
  high_mito_flag

colData(sce_vcm)$qc_pass <-
  qc_pass


sce_vcm_qc <- sce_vcm[
  ,
  colData(sce_vcm)$qc_pass
]


dim(sce_vcm_qc)

table(
  colData(sce_vcm_qc)$age_group
)

length(
  unique(colData(sce_vcm_qc)$Donor)
)
# ------------------------------------------------------------
# RECOVERY STEP 6: Save post-QC checkpoint
# ------------------------------------------------------------
#
# This saves the QC-filtered SCE BEFORE normalization/PCA.
#
# If R crashes later, we can reload this object instead of
# repeating the QC workflow.

saveRDS(
  sce_vcm_qc,
  file = "data/sce_vcm_qc_rawcounts.rds",
  compress = FALSE
)

file.exists(
  "data/sce_vcm_qc_rawcounts.rds"
)

file.info(
  "data/sce_vcm_qc_rawcounts.rds"
)$size / 1024^3

# ------------------------------------------------------------
# RECOVERY STEP 7: Calculate library-size factors
# ------------------------------------------------------------

sce_vcm_qc <- scuttle::computeLibraryFactors(
  sce_vcm_qc
)


summary(
  sizeFactors(sce_vcm_qc)
)
# ------------------------------------------------------------
# RECOVERY STEP 8: Log-normalize expression
# ------------------------------------------------------------
#
# Raw counts remain in:
#   assay(sce_vcm_qc, "counts")
#
# Normalized log-expression will be stored separately in:
#   assay(sce_vcm_qc, "logcounts")
#
# We will NOT use logcounts for pseudobulk DE later.

sce_vcm_qc <- scuttle::logNormCounts(
  sce_vcm_qc
)


assayNames(sce_vcm_qc)

# ------------------------------------------------------------
# RECOVERY STEP 9: Model gene variance
# ------------------------------------------------------------

library(scran)


gene_variance <- scran::modelGeneVar(
  sce_vcm_qc
)

summary(
  gene_variance$bio
)
# ------------------------------------------------------------
# RECOVERY STEP 10: Select the top 2,000 HVGs
# ------------------------------------------------------------

n_hvg <- min(
  2000,
  nrow(gene_variance)
)


hvg_order <- order(
  gene_variance$bio,
  decreasing = TRUE
)


hvg_genes <- rownames(
  gene_variance
)[
  hvg_order[
    seq_len(n_hvg)
  ]
]


length(hvg_genes)

head(
  hvg_genes,
  20
)
# ------------------------------------------------------------
# RECOVERY STEP 11: Save normalized pre-PCA checkpoint
# ------------------------------------------------------------
#
# This checkpoint contains:
#
# - 70,364 QC-passing nuclei
# - raw counts
# - log-normalized counts
# - QC metadata
#
# We separately save the HVG list and variance model.
#
# Therefore, if PCA fails again, we can reload everything
# directly instead of repeating QC and normalization.

saveRDS(
  sce_vcm_qc,
  file = "data/sce_vcm_qc_prePCA.rds",
  compress = FALSE
)


saveRDS(
  hvg_genes,
  file = "data/hvg_genes.rds",
  compress = FALSE
)


saveRDS(
  gene_variance,
  file = "data/gene_variance.rds",
  compress = FALSE
)

file.exists(
  c(
    "data/sce_vcm_qc_prePCA.rds",
    "data/hvg_genes.rds",
    "data/gene_variance.rds"
  )
)
sce_vcm_qc <- readRDS(
  "data/sce_vcm_qc_prePCA.rds"
)

hvg_genes <- readRDS(
  "data/hvg_genes.rds"
)

gene_variance <- readRDS(
  "data/gene_variance.rds"
)
# ------------------------------------------------------------
# RECOVERY STEP 12: Prepare HVG matrix for PCA
# ------------------------------------------------------------

hvg_logcounts <- logcounts(sce_vcm_qc)[
  hvg_genes,
  ,
  drop = FALSE
]


cat(
  "Dimensions:",
  dim(hvg_logcounts),
  "\n"
)

cat(
  "Class:",
  class(hvg_logcounts),
  "\n"
)

cat(
  "Memory:",
  format(
    object.size(hvg_logcounts),
    units = "MB"
  ),
  "\n"
)

# ------------------------------------------------------------
# STEP 56: Convert HVG matrix to compressed sparse format
# ------------------------------------------------------------
#
# Current class:
#   dgTMatrix  = sparse triplet matrix
#
# Convert to:
#   dgCMatrix  = compressed sparse column matrix
#
# This does NOT densify the matrix.
# It remains sparse and is better suited for irlba.

hvg_logcounts_csc <- methods::as(
  hvg_logcounts,
  "CsparseMatrix"
)


# Confirm

dim(hvg_logcounts_csc)

class(hvg_logcounts_csc)

format(
  object.size(hvg_logcounts_csc),
  units = "MB"
)

# ------------------------------------------------------------
# STEP 57: Run truncated PCA with irlba
# ------------------------------------------------------------
#
# Input after transpose:
#   rows    = nuclei
#   columns = HVGs
#
# We calculate only the first 30 PCs.
#
# IMPORTANT:
# This should avoid the sparse-to-dense conversion that
# happened with rsvd::rpca().

library(irlba)

set.seed(12345)

pca_result <- irlba::prcomp_irlba(
  t(hvg_logcounts_csc),
  n = 30,
  center = TRUE,
  scale. = FALSE
)
dim(pca_result$x)

head(pca_result$sdev)

# ------------------------------------------------------------
# STEP 58: Store PCA coordinates
# ------------------------------------------------------------

reducedDim(
  sce_vcm_qc,
  "PCA"
) <- pca_result$x


dim(
  reducedDim(
    sce_vcm_qc,
    "PCA"
  )
)
# ------------------------------------------------------------
# STEP 59: Save post-PCA checkpoint
# ------------------------------------------------------------

saveRDS(
  sce_vcm_qc,
  "data/sce_vcm_qc_PCA.rds",
  compress = FALSE
)

saveRDS(
  pca_result,
  "data/pca_result.rds",
  compress = FALSE
)
dim(pca_result$x)
# ------------------------------------------------------------
# STEP 58: Store PCA coordinates in the SCE
# ------------------------------------------------------------

reducedDim(
  sce_vcm_qc,
  "PCA"
) <- pca_result$x


# Confirm

reducedDimNames(
  sce_vcm_qc
)

dim(
  reducedDim(
    sce_vcm_qc,
    "PCA"
  )
)
# ------------------------------------------------------------
# STEP 59: Save post-PCA checkpoint
# ------------------------------------------------------------
#
# This protects us from having to repeat PCA if R crashes.

saveRDS(
  sce_vcm_qc,
  file = "data/sce_vcm_qc_PCA.rds",
  compress = FALSE
)

saveRDS(
  pca_result,
  file = "data/pca_result.rds",
  compress = FALSE
)
file.exists(
  c(
    "data/sce_vcm_qc_PCA.rds",
    "data/pca_result.rds"
  )
)
# ------------------------------------------------------------
# STEP 60: Create PCA plotting table
# ------------------------------------------------------------

pca_plot_data <- as.data.frame(
  reducedDim(
    sce_vcm_qc,
    "PCA"
  )
)

colnames(pca_plot_data) <- paste0(
  "PC",
  seq_len(ncol(pca_plot_data))
)


pca_plot_data <- dplyr::bind_cols(
  pca_plot_data,
  as.data.frame(colData(sce_vcm_qc))
)
# ------------------------------------------------------------
# STEP 61: PCA colored by source study
# ------------------------------------------------------------

p_pca_study <- ggplot(
  pca_plot_data,
  aes(
    x = PC1,
    y = PC2,
    color = DataSource
  )
) +
  geom_point(
    alpha = 0.35,
    size = 0.4
  ) +
  labs(
    title = "VCM PCA by source study",
    x = "PC1",
    y = "PC2",
    color = "Study"
  ) +
  theme_classic()

p_pca_study
# ------------------------------------------------------------
# STEP 62: PCA colored by age group
# ------------------------------------------------------------

p_pca_age <- ggplot(
  pca_plot_data,
  aes(
    x = PC1,
    y = PC2,
    color = age_group
  )
) +
  geom_point(
    alpha = 0.35,
    size = 0.4
  ) +
  labs(
    title = "VCM PCA by female age group",
    x = "PC1",
    y = "PC2",
    color = "Age group"
  ) +
  theme_classic()

p_pca_age
# ------------------------------------------------------------
# STEP 63: Quantify PC association with study and age
# ------------------------------------------------------------
#
# These R-squared values are exploratory.
# They tell us whether each PC is more strongly associated
# with study/batch or with age group.

pc_association_table <- purrr::map_dfr(
  paste0("PC", 1:10),
  function(pc_name) {
    
    study_model <- stats::lm(
      pca_plot_data[[pc_name]] ~
        pca_plot_data$DataSource
    )
    
    age_model <- stats::lm(
      pca_plot_data[[pc_name]] ~
        pca_plot_data$age_group
    )
    
    tibble::tibble(
      PC = pc_name,
      study_r2 = summary(study_model)$r.squared,
      age_group_r2 = summary(age_model)$r.squared
    )
  }
)

print(
  pc_association_table,
  n = Inf
)  
# ------------------------------------------------------------
# STEP 64: Inspect PCA by source study
# ------------------------------------------------------------

p_pca_study

# ------------------------------------------------------------
# STEP 65: PCA by study, colored by age group
# ------------------------------------------------------------
#
# Each panel represents one source study.
#
# This is especially important because age group and study
# are partially confounded in our 24-donor cohort.
#
# We want to ask:
#
# Within studies that contain BOTH younger and older female
# donors, do the age groups still show some separation?
#
# This is much more informative than looking at the pooled
# younger-vs-older PCA alone.

p_pca_age_by_study <- ggplot(
  pca_plot_data,
  aes(
    x = PC1,
    y = PC2,
    color = age_group
  )
) +
  geom_point(
    alpha = 0.35,
    size = 0.35
  ) +
  facet_wrap(
    ~ DataSource,
    scales = "free"
  ) +
  labs(
    title = "VCM PCA by age group within source study",
    subtitle = "Unintegrated PCA",
    x = "PC1",
    y = "PC2",
    color = "Age group"
  ) +
  theme_classic()

p_pca_age_by_study
# ------------------------------------------------------------
# STEP 66: Inspect PCs most associated with age group
# ------------------------------------------------------------

p_pca_pc3_pc7 <- ggplot(
  pca_plot_data,
  aes(
    x = PC3,
    y = PC7,
    color = age_group
  )
) +
  geom_point(
    alpha = 0.35,
    size = 0.35
  ) +
  labs(
    title = "VCM PCA: PC3 versus PC7",
    subtitle = "PCs showing the strongest age-group associations",
    x = "PC3",
    y = "PC7",
    color = "Age group"
  ) +
  theme_classic()

p_pca_pc3_pc7
# ------------------------------------------------------------
# STEP 67: PC3 vs PC7 within each study
# ------------------------------------------------------------

p_pca_pc3_pc7_study <- ggplot(
  pca_plot_data,
  aes(
    x = PC3,
    y = PC7,
    color = age_group
  )
) +
  geom_point(
    alpha = 0.35,
    size = 0.35
  ) +
  facet_wrap(
    ~ DataSource,
    scales = "free"
  ) +
  labs(
    title = "VCM PC3 versus PC7 by source study",
    subtitle = "Exploratory assessment of age-associated structure",
    x = "PC3",
    y = "PC7",
    color = "Age group"
  ) +
  theme_classic()

p_pca_pc3_pc7_study
# ------------------------------------------------------------
# STEP 68: Prepare for UMAP
# ------------------------------------------------------------
#
# We deliberately run UMAP on the UNINTEGRATED PCA.
#
# Purpose:
#   - visualize transcriptional structure
#   - examine study effects
#   - examine donor effects
#   - examine age-group distribution
#
# UMAP is NOT used for differential expression.
#
# We use the existing 30-PC representation rather than
# returning to the expression matrix.

if (!requireNamespace("uwot", quietly = TRUE)) {
# install.packages("uwot")
}

library(uwot)

# ------------------------------------------------------------
# STEP 69: Extract PCA coordinates for UMAP
# ------------------------------------------------------------

pca_for_umap <- reducedDim(
  sce_vcm_qc,
  "PCA"
)

dim(pca_for_umap)

# ------------------------------------------------------------
# STEP 70: Run unintegrated UMAP
# ------------------------------------------------------------
#
# We use all 70,364 QC-passing VCM nuclei.
#
# Input:
#   70,364 nuclei
#   30 PCs
#
# n_neighbors = 30:
#   standard starting point for local neighborhood structure
#
# min_dist = 0.3:
#   moderate compression of local neighborhoods
#
# metric = "cosine":
#   commonly useful for reduced-dimensional single-cell data
#
# n_threads = 1:
#   keeps the result reproducible with set.seed().
#
# We are intentionally NOT correcting for study here.

set.seed(12345)

umap_result <- uwot::umap(
  pca_for_umap,
  n_neighbors = 30,
  min_dist = 0.3,
  metric = "cosine",
  n_components = 2,
  n_threads = 1,
  verbose = TRUE
)
# ------------------------------------------------------------
# STEP 71: Check and store UMAP
# ------------------------------------------------------------

dim(umap_result)

head(umap_result)
reducedDim(
  sce_vcm_qc,
  "UMAP_unintegrated"
) <- umap_result


dim(
  reducedDim(
    sce_vcm_qc,
    "UMAP_unintegrated"
  )
)

# ------------------------------------------------------------
# STEP 72: Save UMAP checkpoint
# ------------------------------------------------------------

saveRDS(
  sce_vcm_qc,
  file = "data/sce_vcm_qc_PCA_UMAP.rds",
  compress = FALSE
)

saveRDS(
  umap_result,
  file = "data/umap_unintegrated.rds",
  compress = FALSE
)
# ------------------------------------------------------------
# STEP 73: Build UMAP plotting table
# ------------------------------------------------------------

umap_plot_data <- as.data.frame(
  reducedDim(
    sce_vcm_qc,
    "UMAP_unintegrated"
  )
)

colnames(umap_plot_data) <- c(
  "UMAP1",
  "UMAP2"
)


umap_plot_data <- dplyr::bind_cols(
  umap_plot_data,
  as.data.frame(colData(sce_vcm_qc))
)
# ------------------------------------------------------------
# STEP 74: UMAP by source study
# ------------------------------------------------------------

p_umap_study <- ggplot(
  umap_plot_data,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = DataSource
  )
) +
  geom_point(
    size = 0.35,
    alpha = 0.35
  ) +
  labs(
    title = "Unintegrated VCM UMAP by source study",
    color = "Study"
  ) +
  theme_classic()

p_umap_study

# ------------------------------------------------------------
# STEP 75: UMAP by female age group
# ------------------------------------------------------------

p_umap_age <- ggplot(
  umap_plot_data,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = age_group
  )
) +
  geom_point(
    size = 0.35,
    alpha = 0.35
  ) +
  labs(
    title = "Unintegrated VCM UMAP by female age group",
    color = "Age group"
  ) +
  theme_classic()

p_umap_age

# ------------------------------------------------------------
# STEP 76: Prepare Harmony batch correction
# ------------------------------------------------------------
#
# WHY:
#
# The unintegrated PCA and UMAP show strong source-study
# structure.
#
# We therefore create a batch-corrected representation
# specifically for:
#
#   - visualization
#   - neighborhood construction
#   - possible VCM subclustering
#
# IMPORTANT:
#
# Harmony-corrected coordinates will NOT be used for
# pseudobulk differential expression.
#
# We correct:
#   DataSource
#
# We DO NOT correct:
#   age_group
#
# because age is the biological variable of interest.

if (!requireNamespace("harmony", quietly = TRUE)) {
# install.packages("harmony")
}

library(harmony)

# ------------------------------------------------------------
# STEP 77: Confirm study × age structure
# ------------------------------------------------------------

study_age_before_harmony <- table(
  colData(sce_vcm_qc)$DataSource,
  colData(sce_vcm_qc)$age_group
)

study_age_before_harmony

# ------------------------------------------------------------
# ------------------------------------------------------------
# STEP 78: Harmony correction of PCA coordinates
# ------------------------------------------------------------
#
# Your installed Harmony version uses RunHarmony()
# rather than the older HarmonyMatrix() function.
#
# Input:
#   70,364 nuclei x 30 PCs
#
# Batch variable:
#   DataSource
#
# IMPORTANT:
#
# We correct ONLY the reduced-dimensional PCA representation.
#
# We DO NOT modify:
#   - raw counts
#   - logcounts
#   - age_group
#
# Harmony will be used only for visualization / possible
# subclustering, NOT for pseudobulk differential expression.
# ------------------------------------------------------------

set.seed(12345)

harmony_result <- harmony::RunHarmony(
  
  data_mat = reducedDim(
    sce_vcm_qc,
    "PCA"
  ),
  
  meta_data = as.data.frame(
    colData(sce_vcm_qc)
  ),
  
  vars_use = "DataSource",
  
  verbose = TRUE
)

# ------------------------------------------------------------
# STEP 79: Inspect Harmony result
# ------------------------------------------------------------

dim(harmony_result)

class(harmony_result)

# ------------------------------------------------------------
# STEP 80: Store Harmony coordinates in the SCE
# ------------------------------------------------------------

reducedDim(
  sce_vcm_qc,
  "HARMONY"
) <- harmony_result


reducedDimNames(
  sce_vcm_qc
)


dim(
  reducedDim(
    sce_vcm_qc,
    "HARMONY"
  )
)

# ------------------------------------------------------------
# STEP 81: Save post-Harmony checkpoint
# ------------------------------------------------------------

saveRDS(
  sce_vcm_qc,
  file = "data/sce_vcm_qc_PCA_UMAP_Harmony.rds",
  compress = FALSE
)


saveRDS(
  harmony_result,
  file = "data/harmony_result.rds",
  compress = FALSE
)

# ------------------------------------------------------------
# STEP 82: Run UMAP using Harmony-corrected coordinates
# ------------------------------------------------------------
#
# Input:
#   70,364 nuclei x 30 Harmony dimensions
#
# This UMAP is for visualization only.
#
# We are NOT using Harmony values for pseudobulk DE.

set.seed(12345)

umap_harmony <- uwot::umap(
  reducedDim(
    sce_vcm_qc,
    "HARMONY"
  ),
  n_neighbors = 30,
  min_dist = 0.3,
  metric = "cosine",
  n_components = 2,
  n_threads = 1,
  verbose = TRUE
)
# ------------------------------------------------------------
# STEP 83: Store Harmony UMAP
# ------------------------------------------------------------

dim(umap_harmony)

reducedDim(
  sce_vcm_qc,
  "UMAP_HARMONY"
) <- umap_harmony


reducedDimNames(
  sce_vcm_qc
)

# ------------------------------------------------------------
# STEP 84: Save Harmony UMAP checkpoint
# ------------------------------------------------------------

saveRDS(
  sce_vcm_qc,
  file = "data/sce_vcm_qc_PCA_UMAP_Harmony_final.rds",
  compress = FALSE
)

saveRDS(
  umap_harmony,
  file = "data/umap_harmony.rds",
  compress = FALSE
)
# ------------------------------------------------------------
# STEP 85: Create Harmony UMAP plotting table
# ------------------------------------------------------------

umap_harmony_data <- as.data.frame(
  reducedDim(
    sce_vcm_qc,
    "UMAP_HARMONY"
  )
)

colnames(umap_harmony_data) <- c(
  "UMAP1",
  "UMAP2"
)


umap_harmony_data <- dplyr::bind_cols(
  umap_harmony_data,
  as.data.frame(colData(sce_vcm_qc))
)
# ------------------------------------------------------------
# STEP 86: Harmony UMAP by source study
# ------------------------------------------------------------

p_harmony_study <- ggplot(
  umap_harmony_data,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = DataSource
  )
) +
  geom_point(
    size = 0.35,
    alpha = 0.35
  ) +
  labs(
    title = "Harmony-corrected VCM UMAP by source study",
    color = "Study"
  ) +
  theme_classic()

p_harmony_study
# ------------------------------------------------------------
# STEP 87: Harmony UMAP by female age group
# ------------------------------------------------------------

p_harmony_age <- ggplot(
  umap_harmony_data,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = age_group
  )
) +
  geom_point(
    size = 0.35,
    alpha = 0.35
  ) +
  labs(
    title = "Harmony-corrected VCM UMAP by female age group",
    color = "Age group"
  ) +
  theme_classic()

p_harmony_age

# ------------------------------------------------------------
# STEP 88: Harmony UMAP by age group within study
# ------------------------------------------------------------

p_harmony_age_by_study <- ggplot(
  umap_harmony_data,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = age_group
  )
) +
  geom_point(
    size = 0.30,
    alpha = 0.35
  ) +
  facet_wrap(
    ~ DataSource
  ) +
  labs(
    title = "Harmony-corrected VCM UMAP by age group within study",
    color = "Age group"
  ) +
  theme_classic()

p_harmony_age_by_study

# ------------------------------------------------------------
# STEP 89: Build nearest-neighbor graph from Harmony space
# ------------------------------------------------------------
#
# We construct the graph using the batch-corrected Harmony
# representation.
#
# IMPORTANT:
#
# Harmony is used here ONLY for:
#   - neighborhood construction
#   - clustering
#   - visualization
#
# It will NOT be used for differential expression.
#
# We use the first 30 Harmony dimensions.

library(scran)

set.seed(12345)

vcm_graph <- scran::buildSNNGraph(
  sce_vcm_qc,
  use.dimred = "HARMONY",
  k = 30
)
vcm_graph
# ------------------------------------------------------------
# STEP 90: Cluster the VCM graph
# ------------------------------------------------------------
#
# Leiden community detection identifies groups of nuclei
# occupying similar transcriptional neighborhoods.
#
# These are initially UNSUPERVISED clusters.
#
# We should NOT assign biological names yet.

set.seed(12345)

vcm_cluster <- igraph::cluster_leiden(
  vcm_graph,
  objective_function = "modularity"
)

vcm_cluster
# ------------------------------------------------------------
# STEP 91: Add cluster membership to the SCE
# ------------------------------------------------------------

colData(sce_vcm_qc)$VCM_cluster <- factor(
  igraph::membership(vcm_cluster)
)


table(
  colData(sce_vcm_qc)$VCM_cluster
)
# ------------------------------------------------------------
# STEP 92: Summarize cluster sizes
# ------------------------------------------------------------

cluster_size_table <- as.data.frame(
  colData(sce_vcm_qc)
) %>%
  
  dplyr::count(
    VCM_cluster,
    name = "n_nuclei"
  ) %>%
  
  dplyr::mutate(
    percent =
      100 * n_nuclei / sum(n_nuclei)
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(n_nuclei)
  )


cluster_size_table
# ------------------------------------------------------------
# STEP 93: Harmony UMAP colored by unsupervised VCM cluster
# ------------------------------------------------------------

umap_harmony_data$VCM_cluster <-
  colData(sce_vcm_qc)$VCM_cluster


p_harmony_cluster <- ggplot(
  umap_harmony_data,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = VCM_cluster
  )
) +
  geom_point(
    size = 0.35,
    alpha = 0.45
  ) +
  labs(
    title = "Harmony-corrected VCM UMAP by unsupervised cluster",
    color = "Cluster"
  ) +
  theme_classic()

p_harmony_cluster
# ------------------------------------------------------------
# STEP 94: Cluster composition by source study
# ------------------------------------------------------------

cluster_study_table <- as.data.frame(
  colData(sce_vcm_qc)
) %>%
  
  dplyr::count(
    VCM_cluster,
    DataSource,
    name = "n_nuclei"
  ) %>%
  
  dplyr::group_by(
    VCM_cluster
  ) %>%
  
  dplyr::mutate(
    percent_within_cluster =
      100 * n_nuclei / sum(n_nuclei)
  ) %>%
  
  dplyr::ungroup()


print(
  cluster_study_table,
  n = Inf
)
# ------------------------------------------------------------
# STEP 95: Cluster composition by age group
# ------------------------------------------------------------
#
# IMPORTANT:
#
# This is descriptive only.
#
# We cannot use nucleus counts as independent biological
# replicates to claim an age-related composition difference.

cluster_age_table <- as.data.frame(
  colData(sce_vcm_qc)
) %>%
  
  dplyr::count(
    VCM_cluster,
    age_group,
    name = "n_nuclei"
  ) %>%
  
  dplyr::group_by(
    VCM_cluster
  ) %>%
  
  dplyr::mutate(
    percent_within_cluster =
      100 * n_nuclei / sum(n_nuclei)
  ) %>%
  
  dplyr::ungroup()


print(
  cluster_age_table,
  n = Inf
)
# ------------------------------------------------------------
# STEP 96: Define canonical ventricular cardiomyocyte markers
# ------------------------------------------------------------
#
# The source datasets already supplied harmonized cell-type
# annotations and our analysis is restricted to cells labeled
# Ventricular_Cardiomyocytes.
#
# We therefore do NOT perform de novo whole-heart annotation.
#
# Instead, we verify that the retained nuclei show the expected
# cardiomyocyte transcriptional identity.

vcm_marker_genes <- c(
  "TTN",
  "MYH7",
  "MYH6",
  "TNNT2",
  "TNNI3",
  "ACTC1",
  "MYL2",
  "MYL3",
  "RYR2",
  "PLN",
  "ATP2A2",
  "CASQ2"
)


# Check which markers are available in this shared-gene matrix

vcm_marker_check <- tibble::tibble(
  gene = vcm_marker_genes,
  present =
    vcm_marker_genes %in%
    rowData(sce_vcm_qc)$gene_short_name
)

vcm_marker_check
# ------------------------------------------------------------
# STEP 97: Summarize canonical VCM marker expression
# ------------------------------------------------------------
#
# For each marker calculate:
#
#   1. total raw counts
#   2. number of nuclei with detectable expression
#   3. percentage of nuclei with detectable expression
#   4. mean normalized log-expression
#
# These values provide a simple identity check.

marker_rows <- match(
  vcm_marker_genes,
  rowData(sce_vcm_qc)$gene_short_name
)


# Keep only genes actually present in the dataset

marker_rows <- marker_rows[
  !is.na(marker_rows)
]


vcm_marker_summary <- tibble::tibble(
  
  gene =
    rowData(sce_vcm_qc)$gene_short_name[
      marker_rows
    ],
  
  total_counts =
    Matrix::rowSums(
      counts(sce_vcm_qc)[
        marker_rows,
        ,
        drop = FALSE
      ]
    ),
  
  nuclei_detected =
    Matrix::rowSums(
      counts(sce_vcm_qc)[
        marker_rows,
        ,
        drop = FALSE
      ] > 0
    ),
  
  percent_nuclei_detected =
    100 *
    nuclei_detected /
    ncol(sce_vcm_qc),
  
  mean_log_expression =
    Matrix::rowMeans(
      logcounts(sce_vcm_qc)[
        marker_rows,
        ,
        drop = FALSE
      ]
    )
)


print(
  vcm_marker_summary,
  n = Inf
)
# ------------------------------------------------------------
# STEP 98: Define broad non-cardiomyocyte markers
# ------------------------------------------------------------
#
# These are NOT being used to re-annotate every nucleus.
#
# They are a contamination sanity check for major:
#
#   fibroblast
#   endothelial
#   immune
#   perivascular
#   adipocyte
#
# signatures.

non_vcm_marker_genes <- c(
  
  # Fibroblast
  "COL1A1",
  "COL1A2",
  "DCN",
  "LUM",
  
  # Endothelial
  "PECAM1",
  "VWF",
  "KDR",
  "EMCN",
  
  # Immune
  "PTPRC",
  "CD68",
  "LYZ",
  
  # Perivascular / smooth muscle
  "RGS5",
  "PDGFRB",
  "ACTA2",
  
  # Adipocyte
  "ADIPOQ",
  "PLIN1"
)


non_vcm_marker_check <- tibble::tibble(
  
  gene = non_vcm_marker_genes,
  
  present =
    non_vcm_marker_genes %in%
    rowData(sce_vcm_qc)$gene_short_name
)


non_vcm_marker_check

# ------------------------------------------------------------
# STEP 99: Save exploratory cluster summaries
# ------------------------------------------------------------

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)


readr::write_csv(
  cluster_size_table,
  "results/tables/VCM_exploratory_cluster_sizes.csv"
)


readr::write_csv(
  cluster_study_table,
  "results/tables/VCM_exploratory_clusters_by_study.csv"
)


readr::write_csv(
  cluster_age_table,
  "results/tables/VCM_exploratory_clusters_by_age.csv"
)
# ------------------------------------------------------------
# INTERPRETATION NOTE
# ------------------------------------------------------------
#
# Harmony substantially improved cross-study mixing compared
# with the unintegrated PCA/UMAP.
#
# However, many Leiden clusters remained strongly dominated
# by individual source studies.
#
# Examples include clusters dominated by Chaffin, Koenig,
# Litvinukova, or Tucker.
#
# Therefore these clusters are retained as EXPLORATORY
# transcriptional neighborhoods and are NOT interpreted as
# validated biological VCM subtypes.
#
# Apparent age enrichment of individual clusters is also not
# interpreted directly because age group and source study are
# partially confounded.
#
# Primary age-associated transcriptional inference will use
# donor-level pseudobulk counts with study-aware modeling.

