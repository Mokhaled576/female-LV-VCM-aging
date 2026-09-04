# ============================================================
# 04_VCM_pseudobulk_age_DE.R
# ============================================================
#
# PROJECT:
# Female human ventricular cardiomyocyte aging analysis
#
# PRIMARY QUESTION:
#
# What transcriptional programs distinguish ventricular
# cardiomyocytes from younger (20-45 years) and older
# (>=55 years) female donors in the non-failing human heart?
#
# Particular biological interests:
#
#   - CYP-mediated metabolism
#   - EET / HETE / arachidonic acid pathways
#   - soluble epoxide hydrolase (EPHX2 / sEH)
#   - estrogen-associated signaling
#   - PPAR / PGC-1 metabolic signaling
#   - AHR signaling
#   - oxidative / metabolic pathways
#
#
# IMPORTANT STATISTICAL PRINCIPLE:
#
# The biological replicate is the DONOR, not the nucleus.
#
# Therefore:
#
#   70,364 nuclei != 70,364 independent replicates
#
# Instead, raw counts from nuclei belonging to the same donor
# will be summed to produce one pseudobulk library per donor.
#
# Final design:
#
#   24 female donors
#
#     5 younger: 20-45 years
#    19 older:   >=55 years
#
# All nuclei:
#
#   ventricular cardiomyocytes
#   left ventricle
#   QC-passing
#
# ============================================================
# ------------------------------------------------------------
# STEP 1: Clean environment
# ------------------------------------------------------------
#
# Script 04 should be reproducible independently of Script 03.
#
# All required information will be loaded from the saved
# post-QC/preprocessing checkpoint.

rm(list = ls())

gc()

set.seed(12345)

# ------------------------------------------------------------
# STEP 2: Load required packages
# ------------------------------------------------------------

library(SingleCellExperiment)
library(Matrix)
library(edgeR)

# Load tidyverse LAST to reduce function masking problems.
library(tidyverse)

# ------------------------------------------------------------
# STEP 3: Load QC-passing VCM object
# ------------------------------------------------------------

sce_file <-
  "data/sce_vcm_qc_PCA_UMAP_Harmony_final.rds"


stopifnot(
  file.exists(sce_file)
)


sce_vcm_qc <- readRDS(
  sce_file
)

# ------------------------------------------------------------
# STEP 4: Verify the loaded dataset
# ------------------------------------------------------------

cat(
  "Genes:",
  nrow(sce_vcm_qc),
  "\n"
)

cat(
  "QC-passing nuclei:",
  ncol(sce_vcm_qc),
  "\n"
)

cat(
  "Donors:",
  length(unique(colData(sce_vcm_qc)$Donor)),
  "\n"
)

cat(
  "Studies:",
  length(unique(colData(sce_vcm_qc)$DataSource)),
  "\n"
)


table(
  colData(sce_vcm_qc)$age_group
)

counts(sce_vcm_qc)
logcounts(sce_vcm_qc)

# ------------------------------------------------------------
# STEP 5: Verify raw-count assay
# ------------------------------------------------------------

assayNames(
  sce_vcm_qc
)


raw_counts <- counts(
  sce_vcm_qc
)


dim(
  raw_counts
)

class(
  raw_counts
)
# ------------------------------------------------------------
# ------------------------------------------------------------
# ------------------------------------------------------------
# STEP 6: Count QC-passing VCM nuclei per donor
# ------------------------------------------------------------
#
# We summarize the number of QC-passing ventricular
# cardiomyocyte nuclei contributed by each donor.
#
# This is descriptive only.
# Donors, not nuclei, will be the biological replicates
# for pseudobulk differential expression.

donor_nuclei_table <- as.data.frame(
  colData(sce_vcm_qc)
) %>%
  
  dplyr::count(
    DataSource,
    Donor,
    Age,
    Sex,
    age_group,
    name = "n_qc_vcm_nuclei"
  ) %>%
  
  dplyr::arrange(
    age_group,
    Age
  )


# Display all 24 donors
donor_nuclei_table

# ------------------------------------------------------------
# STEP 6A: Verify donor summary
# ------------------------------------------------------------

cat(
  "Number of donors:",
  nrow(donor_nuclei_table),
  "\n"
)

cat(
  "Total QC-passing nuclei:",
  sum(donor_nuclei_table$n_qc_vcm_nuclei),
  "\n"
)
# ------------------------------------------------------------
# STEP 7: Define donor membership for every nucleus
# ------------------------------------------------------------
#
# Each nucleus is assigned to one of the 24 donors.
#
# factor() gives us a consistent donor ordering that will
# also be used for the pseudobulk count matrix.

donor_factor <- factor(
  colData(sce_vcm_qc)$Donor
)


levels(donor_factor)

length(levels(donor_factor))

# ------------------------------------------------------------
# STEP 8: Build sparse nucleus-to-donor aggregation matrix
# ------------------------------------------------------------
#
# Rows:
#   70,364 nuclei
#
# Columns:
#   24 donors
#
# Each row contains a 1 in the column representing the
# donor from whom that nucleus originated.
#
# The matrix stays sparse and is therefore memory efficient.

donor_design <- Matrix::sparse.model.matrix(
  ~ 0 + donor_factor
)


# Rename columns using the actual donor IDs

colnames(donor_design) <- levels(
  donor_factor
)


dim(donor_design)

# ------------------------------------------------------------
# STEP 9: Aggregate raw nucleus counts by donor
# ------------------------------------------------------------
#
# IMPORTANT:
#
# We use RAW counts:
#
#   raw_counts = genes x nuclei
#
# and multiply by:
#
#   donor_design = nuclei x donors
#
# producing:
#
#   pb_counts = genes x donors
#
# Thus every donor becomes one independent pseudobulk
# RNA-seq-like library.

pb_counts <- raw_counts %*% donor_design


dim(pb_counts)

class(pb_counts)

# ------------------------------------------------------------
# STEP 10: Verify conservation of total raw counts
# ------------------------------------------------------------
#
# Pseudobulk aggregation should only SUM counts.
#
# Therefore the total number of counts before and after
# aggregation must be identical.

total_cell_counts <- sum(
  raw_counts
)


total_pseudobulk_counts <- sum(
  pb_counts
)


cat(
  "Total counts across nuclei:",
  total_cell_counts,
  "\n"
)

cat(
  "Total counts after pseudobulk:",
  total_pseudobulk_counts,
  "\n"
)


all.equal(
  total_cell_counts,
  total_pseudobulk_counts
)

# ------------------------------------------------------------
# STEP 11: Construct donor-level pseudobulk metadata
# ------------------------------------------------------------
#
# The count matrix now contains one column per donor,
# therefore our metadata must also contain exactly one
# row per donor.

pb_meta <- as.data.frame(
  colData(sce_vcm_qc)
) %>%
  
  dplyr::select(
    Donor,
    DataSource,
    Age,
    Sex,
    age_group
  ) %>%
  
  dplyr::distinct()


# Check number of unique donors

nrow(pb_meta)
# ------------------------------------------------------------
# STEP 12: Align donor metadata to pseudobulk columns
# ------------------------------------------------------------
#
# This is critical.
#
# Row 1 of pb_meta must correspond to column 1 of pb_counts,
# row 2 to column 2, and so on.

pb_meta <- pb_meta[
  match(
    colnames(pb_counts),
    pb_meta$Donor
  ),
  ,
  drop = FALSE
]


rownames(pb_meta) <- pb_meta$Donor

identical(
  rownames(pb_meta),
  colnames(pb_counts)
)

sum(
  is.na(pb_meta$Donor)
)
# ------------------------------------------------------------
# STEP 13: Add QC-passing nucleus counts to donor metadata
# ------------------------------------------------------------

nuclei_per_donor <- table(
  colData(sce_vcm_qc)$Donor
)


pb_meta$n_qc_vcm_nuclei <- as.integer(
  nuclei_per_donor[
    pb_meta$Donor
  ]
)
sum(
  pb_meta$n_qc_vcm_nuclei
)

# ------------------------------------------------------------
# STEP 14: Calculate pseudobulk library metrics
# ------------------------------------------------------------
#
# library_size:
#   total raw counts in each donor pseudobulk
#
# detected_genes:
#   number of genes having at least one count in that donor

pb_meta$library_size <- Matrix::colSums(
  pb_counts
)


pb_meta$detected_genes <- Matrix::colSums(
  pb_counts > 0
)

# ------------------------------------------------------------
# STEP 14A: Inspect donor pseudobulk libraries
# ------------------------------------------------------------

pb_meta_display <- pb_meta %>%
  
  dplyr::select(
    Donor,
    DataSource,
    Age,
    age_group,
    n_qc_vcm_nuclei,
    library_size,
    detected_genes
  ) %>%
  
  dplyr::arrange(
    library_size
  )


# Convert explicitly to tibble so print(n = Inf) works safely

pb_meta_display <- tibble::as_tibble(
  pb_meta_display
)


print(
  pb_meta_display,
  n = Inf
)
# ------------------------------------------------------------
# STEP 15: Study x age-group donor structure
# ------------------------------------------------------------
#
# This reveals the partial confounding between source study
# and age group.
#
# It will be critical when we construct the edgeR design.

study_age_donor_table <- pb_meta %>%
  
  dplyr::count(
    DataSource,
    age_group,
    name = "n_donors"
  )


study_age_donor_table <- tibble::as_tibble(
  study_age_donor_table
)


print(
  study_age_donor_table,
  n = Inf
)
table(
  pb_meta$DataSource,
  pb_meta$age_group
)
# ------------------------------------------------------------
# STEP 16: Save QC-filtered raw pseudobulk
# ------------------------------------------------------------

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)


saveRDS(
  pb_counts,
  file = "data/female_LV_VCM_QC_pseudobulk_counts.rds",
  compress = FALSE
)


readr::write_csv(
  pb_meta,
  file =
    "results/tables/female_LV_VCM_QC_pseudobulk_metadata.csv"
)


readr::write_csv(
  donor_nuclei_table,
  file =
    "results/tables/female_LV_VCM_QC_nuclei_per_donor.csv"
)

# ------------------------------------------------------------
# STEP 17: Create edgeR DGEList object
# ------------------------------------------------------------
#
# pb_counts contains:
#
#   rows    = genes
#   columns = donor pseudobulk samples
#
# Each donor is now one biological replicate.
#
# IMPORTANT:
# We are using raw summed counts here.
# edgeR requires raw integer-like count data.

dge <- edgeR::DGEList(
  counts = pb_counts,
  samples = pb_meta
)


# Inspect basic structure

dge
# ------------------------------------------------------------
# STEP 18: Verify pseudobulk library sizes
# ------------------------------------------------------------
#
# edgeR calculates library size as the column sum
# of the raw count matrix.
#
# These values should match the library_size column
# we calculated manually.

all.equal(
  as.numeric(dge$samples$lib.size),
  as.numeric(pb_meta$library_size)
)

# ------------------------------------------------------------
# STEP 19: Define age group for edgeR
# ------------------------------------------------------------
#
# Set younger as the reference level.
#
# Therefore later:
#
#   positive logFC = higher expression in older females
#   negative logFC = higher expression in younger females

dge$samples$age_group <- factor(
  dge$samples$age_group,
  levels = c(
    "Younger_20_45",
    "Older_55plus"
  )
)
# ------------------------------------------------------------
# STEP 20: Define source study as a factor
# ------------------------------------------------------------

dge$samples$DataSource <- factor(
  dge$samples$DataSource
)


levels(
  dge$samples$DataSource
)
# ------------------------------------------------------------
# STEP 21: Construct preliminary study + age design
# ------------------------------------------------------------
#
# This design reflects the structure we ultimately want
# to account for:
#
#   source study
#   age group
#
# We are NOT fitting the DE model yet.
#
# First we use the design to perform expression-aware
# gene filtering.

design_prelim <- model.matrix(
  ~ DataSource + age_group,
  data = dge$samples
)


dim(
  design_prelim
)


colnames(
  design_prelim
)


levels(
  dge$samples$age_group
)
# ------------------------------------------------------------
# STEP 22: Check design-matrix rank
# ------------------------------------------------------------
#
# This is critical because age group and DataSource are
# partially confounded.
#
# A full-rank design means the model coefficients are
# mathematically estimable.
#
# It does NOT mean the confounding problem disappears;
# it only means the model can technically be fitted.

design_rank <- qr(
  design_prelim
)$rank


cat(
  "Number of design columns:",
  ncol(design_prelim),
  "\n"
)

cat(
  "Design rank:",
  design_rank,
  "\n"
)

cat(
  "Full rank:",
  design_rank == ncol(design_prelim),
  "\n"
)
# ------------------------------------------------------------
# STEP 23: Filter low-expression genes
# ------------------------------------------------------------
#
# filterByExpr() determines whether each gene has sufficient
# expression across enough donor-level libraries to support
# reliable differential-expression testing.
#
# This is preferable to using an arbitrary raw-count cutoff.
#
# The filtering is performed at the DONOR level.

keep_gene <- edgeR::filterByExpr(
  dge,
  design = design_prelim
)


table(
  keep_gene
)
# Number of genes before filtering

cat(
  "Genes before filtering:",
  nrow(dge),
  "\n"
)


cat(
  "Genes retained:",
  sum(keep_gene),
  "\n"
)


cat(
  "Genes removed:",
  sum(!keep_gene),
  "\n"
)


cat(
  "Percent retained:",
  round(
    100 * mean(keep_gene),
    2
  ),
  "%\n"
)
# ------------------------------------------------------------
# STEP 24: Apply gene-expression filtering
# ------------------------------------------------------------

dge_filtered <- dge[
  keep_gene,
  ,
  keep.lib.sizes = FALSE
]


dim(
  dge_filtered
)
# ------------------------------------------------------------
# STEP 25: TMM normalization
# ------------------------------------------------------------
#
# calcNormFactors() adjusts for differences in library
# composition between donor pseudobulk samples.
#
# It does NOT make the raw library sizes equal.
#
# TMM normalization factors are used internally by edgeR
# when calculating effective library sizes.

dge_filtered <- edgeR::calcNormFactors(
  dge_filtered,
  method = "TMM"
)


dge_filtered$samples %>%
  tibble::as_tibble() %>%
  dplyr::select(
    Donor,
    DataSource,
    Age,
    age_group,
    lib.size,
    norm.factors
  )
# ------------------------------------------------------------
# STEP 26: Calculate effective library sizes
# ------------------------------------------------------------

dge_filtered$samples$effective_lib_size <-
  dge_filtered$samples$lib.size *
  dge_filtered$samples$norm.factors


dge_filtered$samples %>%
  tibble::as_tibble() %>%
  dplyr::select(
    Donor,
    DataSource,
    age_group,
    lib.size,
    norm.factors,
    effective_lib_size
  ) %>%
  dplyr::arrange(
    effective_lib_size
  )
# ------------------------------------------------------------
# STEP 27: Donor-level MDS plot
# ------------------------------------------------------------
#
# This is now much more relevant for our DE question than
# the nucleus-level UMAP.
#
# Each point = one donor.
#
# We first color by age group.


plotMDS(
  dge_filtered,
  labels = dge_filtered$samples$Donor,
  col = as.integer(
    dge_filtered$samples$age_group
  ),
  main = "Donor-level VCM pseudobulk MDS"
)
# ------------------------------------------------------------
# STEP 28: Extract MDS coordinates
# ------------------------------------------------------------

mds_result <- plotMDS(
  dge_filtered,
  plot = FALSE
)



mds_data <- tibble::tibble(
  
  Donor =
    dge_filtered$samples$Donor,
  
  DataSource =
    dge_filtered$samples$DataSource,
  
  Age =
    dge_filtered$samples$Age,
  
  age_group =
    dge_filtered$samples$age_group,
  
  MDS1 =
    mds_result$x,
  
  MDS2 =
    mds_result$y
)
# ------------------------------------------------------------
# STEP 29: Publication-quality donor-level pseudobulk MDS
# ------------------------------------------------------------
#
# Each point represents ONE donor.
#
# Color = age group
# Shape = source study
#
# Donor labels are positioned using ggrepel so that they
# do not overlap as heavily as in the base edgeR plot.


# Install ggrepel only if it is not already installed
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}


library(ggrepel)


p_mds_combined <- ggplot(
  mds_data,
  aes(
    x = MDS1,
    y = MDS2,
    color = age_group,
    shape = DataSource
  )
) +
  
  # Plot donor pseudobulk samples
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  # Add donor IDs while reducing label overlap
  ggrepel::geom_text_repel(
    aes(label = Donor),
    size = 3.5,
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    seed = 12345
  ) +
  
  # Cleaner names for the age groups
  scale_color_discrete(
    labels = c(
      "Younger (20-45)",
      "Older (>=55)"
    )
  ) +
  
  labs(
    title = "Donor-level ventricular cardiomyocyte pseudobulk MDS",
    subtitle = "Female non-failing left ventricle",
    x = "Leading logFC dimension 1",
    y = "Leading logFC dimension 2",
    color = "Age group",
    shape = "Study"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    
    plot.subtitle = element_text(
      size = 12
    ),
    
    legend.position = "right"
  )


p_mds_combined
# ------------------------------------------------------------
# STEP 29B: Refined combined MDS figure
# ------------------------------------------------------------

p_mds_combined <- ggplot(
  mds_data,
  aes(
    x = MDS1,
    y = MDS2,
    color = age_group,
    shape = DataSource
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  ggrepel::geom_text_repel(
    aes(label = Donor),
    size = 3.3,
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 0.7,
    point.padding = 0.4,
    min.segment.length = 0,
    force = 1.5,
    max.time = 2,
    seed = 12345
  ) +
  
  scale_color_discrete(
    labels = c(
      "Younger (20-45)",
      "Older (55+)"
    )
  ) +
  
  labs(
    title = "Donor-level ventricular cardiomyocyte pseudobulk MDS",
    subtitle = "Female non-failing left ventricle",
    x = "Leading logFC dimension 1",
    y = "Leading logFC dimension 2",
    color = "Age group",
    shape = "Study"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      size = 12
    ),
    legend.position = "right"
  )

p_mds_combined
# ------------------------------------------------------------
# STEP 30: MDS colored specifically by source study
# ------------------------------------------------------------
#
# This plot is important because the nucleus-level analysis
# already showed a substantial study effect.
#
# We now check whether that effect is also visible after
# donor-level pseudobulk aggregation.

p_mds_study <- ggplot(
  mds_data,
  aes(
    x = MDS1,
    y = MDS2,
    color = DataSource
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  ggrepel::geom_text_repel(
    aes(label = Donor),
    size = 3.5,
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    seed = 12345
  ) +
  
  labs(
    title = "Donor-level VCM pseudobulk MDS by source study",
    subtitle = "Female non-failing left ventricle",
    x = "Leading logFC dimension 1",
    y = "Leading logFC dimension 2",
    color = "Study"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    
    legend.position = "right"
  )


p_mds_study
# ------------------------------------------------------------
# STEP 30B: Clean study-level pseudobulk MDS
# ------------------------------------------------------------
#
# Each point = one donor pseudobulk sample.
#
# The purpose of this figure is to determine whether source
# study contributes strongly to global transcriptional
# differences between donor samples.
#
# Donor labels are intentionally omitted here because the
# scientific message is the study-level structure.

p_mds_study_clean <- ggplot(
  mds_data,
  aes(
    x = MDS1,
    y = MDS2,
    color = DataSource
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  labs(
    title = "Donor-level VCM pseudobulk MDS by source study",
    subtitle = "Female non-failing left ventricle",
    x = "Leading logFC dimension 1",
    y = "Leading logFC dimension 2",
    color = "Study"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    
    plot.subtitle = element_text(
      size = 12
    ),
    
    legend.position = "right"
  )


p_mds_study_clean

# ------------------------------------------------------------
# STEP 31: MDS colored specifically by age group
# ------------------------------------------------------------

p_mds_age <- ggplot(
  mds_data,
  aes(
    x = MDS1,
    y = MDS2,
    color = age_group
  )
) +
  
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  
  ggrepel::geom_text_repel(
    aes(label = Donor),
    size = 3.5,
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    seed = 12345
  ) +
  
  scale_color_discrete(
    labels = c(
      "Younger (20-45)",
      "Older (>=55)"
    )
  ) +
  
  labs(
    title = "Donor-level VCM pseudobulk MDS by age group",
    subtitle = "Female non-failing left ventricle",
    x = "Leading logFC dimension 1",
    y = "Leading logFC dimension 2",
    color = "Age group"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    
    legend.position = "right"
  )


p_mds_age
# ------------------------------------------------------------
# STEP 32: Construct final study-adjusted design matrix
# ------------------------------------------------------------
#
# PRIMARY MODEL:
#
#   expression ~ DataSource + age_group
#
# DataSource:
#   adjusts for systematic differences among the six
#   source studies.
#
# age_group:
#   estimates the older-versus-younger effect AFTER
#   accounting for source study.
#
# Younger is the reference group, so:
#
#   positive age-group coefficient
#       = higher expression in older females
#
#   negative age-group coefficient
#       = higher expression in younger females


dge_filtered$samples$age_group <- factor(
  dge_filtered$samples$age_group,
  levels = c(
    "Younger_20_45",
    "Older_55plus"
  )
)


dge_filtered$samples$DataSource <- factor(
  dge_filtered$samples$DataSource
)


design <- model.matrix(
  ~ DataSource + age_group,
  data = dge_filtered$samples
)


colnames(design)
dim(design)
# ------------------------------------------------------------
# STEP 33: Formally test design-matrix rank
# ------------------------------------------------------------

design_rank <- qr(design)$rank


cat(
  "Number of donor samples:",
  nrow(design),
  "\n"
)

cat(
  "Number of model coefficients:",
  ncol(design),
  "\n"
)

cat(
  "Design rank:",
  design_rank,
  "\n"
)

cat(
  "Full rank:",
  design_rank == ncol(design),
  "\n"
)
# Show the actual coefficients in the model

colnames(design)
# Confirm study x age structure one final time

table(
  dge_filtered$samples$DataSource,
  dge_filtered$samples$age_group
)
# ------------------------------------------------------------
# STEP 34: Estimate gene-wise biological dispersion
# ------------------------------------------------------------
#
# edgeR models donor-level pseudobulk counts using the
# negative binomial distribution.
#
# Dispersion captures biological variability among donors
# beyond simple Poisson sampling variation.
#
# Our design adjusts for:
#
#   1. source study
#   2. age group
#
# IMPORTANT:
# Each column is one DONOR pseudobulk sample.
# Nuclei are NOT treated as independent replicates.

dge_filtered <- edgeR::estimateDisp(
  dge_filtered,
  design = design,
  robust = TRUE
)
dge_filtered$common.dispersion

summary(
  dge_filtered$tagwise.dispersion
)
# ------------------------------------------------------------
# STEP 35: Inspect biological variability
# ------------------------------------------------------------

plotBCV(
  dge_filtered,
  main = "Biological variation across donor pseudobulk samples"
)
# ------------------------------------------------------------
# STEP 36: Fit edgeR quasi-likelihood model
# ------------------------------------------------------------
#
# The quasi-likelihood framework provides an additional
# layer of uncertainty estimation and is appropriate for
# differential-expression inference with biological
# replicates.
#
# robust = TRUE reduces sensitivity to genes with unusual
# dispersion behavior.

fit <- edgeR::glmQLFit(
  dge_filtered,
  design = design,
  robust = TRUE
)
fit
# ------------------------------------------------------------
# STEP 37: Identify older-versus-younger coefficient
# ------------------------------------------------------------

age_coef <- which(
  colnames(design) == "age_groupOlder_55plus"
)


age_coef

colnames(design)[age_coef]

# ------------------------------------------------------------
# STEP 38: Quasi-likelihood test for age-group effect
# ------------------------------------------------------------
#
# Hypothesis:
#
# H0:
#   no expression difference between older and younger
#   female donors after adjustment for source study
#
# H1:
#   expression differs between older and younger female
#   donors after adjustment for source study.

qlf_age <- edgeR::glmQLFTest(
  fit,
  coef = age_coef
)
qlf_age
# ------------------------------------------------------------
# STEP 39: Extract complete genome-wide DE table
# ------------------------------------------------------------

de_age <- edgeR::topTags(
  qlf_age,
  n = Inf,
  sort.by = "PValue"
)$table


dim(de_age)

head(de_age, 20)
# ------------------------------------------------------------
# STEP 40: Convert DE results into a clean table
# ------------------------------------------------------------

de_age <- de_age %>%
  
  tibble::rownames_to_column(
    var = "gene"
  ) %>%
  
  tibble::as_tibble()


dim(de_age)

# ------------------------------------------------------------
# STEP 41: Summarize genome-wide statistical significance
# ------------------------------------------------------------

de_summary <- de_age %>%
  
  dplyr::summarise(
    
    genes_tested =
      dplyr::n(),
    
    FDR_less_0_05 =
      sum(FDR < 0.05, na.rm = TRUE),
    
    higher_in_older =
      sum(
        FDR < 0.05 &
          logFC > 0,
        na.rm = TRUE
      ),
    
    higher_in_younger =
      sum(
        FDR < 0.05 &
          logFC < 0,
        na.rm = TRUE
      )
  )


de_summary
# ------------------------------------------------------------
# STEP 42: Inspect top age-associated genes
# ------------------------------------------------------------

de_age %>%
  
  dplyr::select(
    gene,
    logFC,
    logCPM,
    F,
    PValue,
    FDR
  ) %>%
  
  dplyr::slice_head(
    n = 30
  )
# ------------------------------------------------------------
# STEP 43: Inspect genome-wide P-value distribution
# ------------------------------------------------------------

p_pvalue_hist <- ggplot(
  de_age,
  aes(x = PValue)
) +
  
  geom_histogram(
    bins = 40,
    boundary = 0
  ) +
  
  labs(
    title = "P-value distribution for the age-group effect",
    subtitle = "Study-adjusted female LV VCM pseudobulk model",
    x = "P-value",
    y = "Number of genes"
  ) +
  
  theme_classic(
    base_size = 13
  )


p_pvalue_hist

# ------------------------------------------------------------
# STEP 44: Save primary all-donor model
# ------------------------------------------------------------

readr::write_csv(
  de_age,
  "results/tables/VCM_age_DE_all24_study_adjusted.csv"
)


readr::write_csv(
  de_summary,
  "results/tables/VCM_age_DE_all24_summary.csv"
)


saveRDS(
  fit,
  "data/VCM_age_edgeR_QL_fit_all24.rds",
  compress = FALSE
)


saveRDS(
  qlf_age,
  "data/VCM_age_edgeR_QL_test_all24.rds",
  compress = FALSE
)
# ------------------------------------------------------------
# STEP 45: Define overlapping-study sensitivity cohort
# ------------------------------------------------------------
#
# PRIMARY ANALYSIS:
#   all 24 donors
#
# SENSITIVITY ANALYSIS:
#   retain only source studies containing BOTH age groups.
#
# These are:
#
#   Koenig
#   Litvinukova
#   Read
#
# This substantially reduces the study-age imbalance.
#
# We intentionally start from the original unfiltered
# donor-level DGEList ("dge"), because expression filtering
# should be recalculated for this smaller cohort.

overlap_studies <- c(
  "Koenig",
  "Litvinukova",
  "Read"
)


keep_overlap_donors <-
  dge$samples$DataSource %in% overlap_studies


table(
  keep_overlap_donors
)
# ------------------------------------------------------------
# STEP 46: Create overlapping-study DGEList
# ------------------------------------------------------------

dge_overlap <- dge[
  ,
  keep_overlap_donors
]


# Drop unused factor levels

dge_overlap$samples$DataSource <- droplevels(
  factor(dge_overlap$samples$DataSource)
)


dge_overlap$samples$age_group <- factor(
  dge_overlap$samples$age_group,
  levels = c(
    "Younger_20_45",
    "Older_55plus"
  )
)
table(
  dge_overlap$samples$DataSource,
  dge_overlap$samples$age_group
)


ncol(
  dge_overlap
)
# ------------------------------------------------------------
# STEP 47: Construct overlapping-study design
# ------------------------------------------------------------
#
# Again:
#
#   expression ~ DataSource + age_group
#
# But now every included study contains both younger
# and older donors.

design_overlap <- model.matrix(
  ~ DataSource + age_group,
  data = dge_overlap$samples
)


colnames(
  design_overlap
)


dim(
  design_overlap
)
# ------------------------------------------------------------
# STEP 48: Check sensitivity-design rank
# ------------------------------------------------------------

cat(
  "Samples:",
  nrow(design_overlap),
  "\n"
)

cat(
  "Coefficients:",
  ncol(design_overlap),
  "\n"
)

cat(
  "Rank:",
  qr(design_overlap)$rank,
  "\n"
)

cat(
  "Full rank:",
  qr(design_overlap)$rank ==
    ncol(design_overlap),
  "\n"
)
# ------------------------------------------------------------
# STEP 49: Recalculate expression filtering
# ------------------------------------------------------------
#
# Do NOT simply inherit the 24-donor filter.
#
# filterByExpr() should reflect the samples included in
# this particular statistical analysis.

keep_gene_overlap <- edgeR::filterByExpr(
  dge_overlap,
  design = design_overlap
)


cat(
  "Genes before filtering:",
  nrow(dge_overlap),
  "\n"
)

cat(
  "Genes retained:",
  sum(keep_gene_overlap),
  "\n"
)

cat(
  "Genes removed:",
  sum(!keep_gene_overlap),
  "\n"
)
dge_overlap <- dge_overlap[
  keep_gene_overlap,
  ,
  keep.lib.sizes = FALSE
]
# ------------------------------------------------------------
# STEP 50: TMM normalization for sensitivity cohort
# ------------------------------------------------------------

dge_overlap <- edgeR::calcNormFactors(
  dge_overlap,
  method = "TMM"
)
# ------------------------------------------------------------
# STEP 51: Fit overlapping-study edgeR model
# ------------------------------------------------------------

dge_overlap <- edgeR::estimateDisp(
  dge_overlap,
  design = design_overlap,
  robust = TRUE
)


fit_overlap <- edgeR::glmQLFit(
  dge_overlap,
  design = design_overlap,
  robust = TRUE
)
# ------------------------------------------------------------
# STEP 52: Test age effect
# ------------------------------------------------------------

age_coef_overlap <- which(
  colnames(design_overlap) ==
    "age_groupOlder_55plus"
)


age_coef_overlap

colnames(design_overlap)[
  age_coef_overlap
]
qlf_overlap <- edgeR::glmQLFTest(
  fit_overlap,
  coef = age_coef_overlap
)
# ------------------------------------------------------------
# STEP 53: Extract sensitivity DE results
# ------------------------------------------------------------

de_overlap <- edgeR::topTags(
  qlf_overlap,
  n = Inf,
  sort.by = "PValue"
)$table %>%
  
  tibble::rownames_to_column(
    var = "gene"
  ) %>%
  
  tibble::as_tibble()
de_overlap_summary <- de_overlap %>%
  
  dplyr::summarise(
    
    genes_tested =
      dplyr::n(),
    
    FDR_less_0_05 =
      sum(
        FDR < 0.05,
        na.rm = TRUE
      ),
    
    higher_in_older =
      sum(
        FDR < 0.05 &
          logFC > 0,
        na.rm = TRUE
      ),
    
    higher_in_younger =
      sum(
        FDR < 0.05 &
          logFC < 0,
        na.rm = TRUE
      )
  )


de_overlap_summary
de_overlap %>%
  
  dplyr::select(
    gene,
    logFC,
    logCPM,
    F,
    PValue,
    FDR
  ) %>%
  
  dplyr::slice_head(
    n = 30
  ) %>%
  
  print(
    n = 30
  )
# ------------------------------------------------------------
# STEP 54: Check CYP4F22 sensitivity
# ------------------------------------------------------------

de_overlap %>%
  
  dplyr::filter(
    gene == "CYP4F22"
  )

# ------------------------------------------------------------
# STEP 55: Compare age-effect estimates between analyses
# ------------------------------------------------------------
#
# We now directly compare logFC estimates from:
#
# PRIMARY:
#   all 24 donors, adjusted for DataSource
#
# SENSITIVITY:
#   15 donors from Koenig, Litvinukova and Read only
#
# The most important question is NOT simply whether each
# analysis independently reaches FDR < 0.05.
#
# We want to know whether estimated age effects have:
#
#   - similar direction
#   - similar magnitude
#
# across the two analysis strategies.


de_comparison <- de_age %>%
  
  dplyr::select(
    gene,
    logFC_primary = logFC,
    PValue_primary = PValue,
    FDR_primary = FDR
  ) %>%
  
  dplyr::inner_join(
    
    de_overlap %>%
      dplyr::select(
        gene,
        logFC_overlap = logFC,
        PValue_overlap = PValue,
        FDR_overlap = FDR
      ),
    
    by = "gene"
  )


dim(de_comparison)
# ------------------------------------------------------------
# STEP 56: Genome-wide effect-size concordance
# ------------------------------------------------------------

cor(
  de_comparison$logFC_primary,
  de_comparison$logFC_overlap,
  method = "pearson"
)


cor(
  de_comparison$logFC_primary,
  de_comparison$logFC_overlap,
  method = "spearman"
)
# ------------------------------------------------------------
# STEP 57: Examine primary FDR-significant genes
#              in the sensitivity analysis
# ------------------------------------------------------------

primary_sig_comparison <- de_comparison %>%
  
  dplyr::filter(
    FDR_primary < 0.05
  )


nrow(
  primary_sig_comparison
)
# ------------------------------------------------------------
# STEP 57B: Account for primary-significant genes that were
#           not tested in the sensitivity analysis
# ------------------------------------------------------------
#
# The primary model identified 111 genes at FDR < 0.05.
#
# However, filterByExpr() was rerun independently for the
# smaller overlapping-study cohort.
#
# Therefore some genes tested in the primary model were
# filtered out and were NOT tested in the sensitivity model.
#
# These genes should be classified as:
#
#   "not tested in sensitivity analysis"
#
# NOT:
#
#   "failed replication"


primary_sig_all <- de_age %>%
  dplyr::filter(
    FDR < 0.05
  )
# ------------------------------------------------------------
# STEP 60: Visualize robustness of primary age effects
# ------------------------------------------------------------
#
# This figure compares the estimated age effect for the
# 95 primary-significant genes that were also testable in
# the overlapping-study sensitivity analysis.
#
# x-axis:
#   log2FC from all 24 donors
#
# y-axis:
#   log2FC from the 15-donor overlapping-study analysis
#
# Points lying close to y = x indicate nearly identical
# effect estimates in the two analyses.


p_logfc_concordance <- ggplot(
  primary_sig_comparison,
  aes(
    x = logFC_primary,
    y = logFC_overlap
  )
) +
  
  # Reference line representing identical effect sizes
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  # Zero-effect reference lines
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  
  geom_point(
    size = 2.8,
    alpha = 0.8
  ) +
  
  labs(
    title = "Age-effect estimates are robust to study restriction",
    subtitle =
      "95 primary-significant genes testable in both analyses",
    x = "log2FC: all 24 donors",
    y = "log2FC: overlapping studies (15 donors)"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )


p_logfc_concordance

nrow(
  primary_sig_all
)
primary_sig_test_status <- primary_sig_all %>%
  
  dplyr::select(
    gene,
    logFC_primary = logFC,
    PValue_primary = PValue,
    FDR_primary = FDR
  ) %>%
  
  dplyr::mutate(
    
    tested_in_overlap =
      gene %in% de_overlap$gene
    
  )
table(
  primary_sig_test_status$tested_in_overlap
)
primary_sig_not_tested <- primary_sig_test_status %>%
  
  dplyr::filter(
    !tested_in_overlap
  ) %>%
  
  dplyr::arrange(
    FDR_primary
  )


print(
  primary_sig_not_tested,
  n = Inf
)
# ------------------------------------------------------------
# STEP 58: Directional concordance among genes tested
#          in BOTH analyses
# ------------------------------------------------------------

primary_sig_comparison <- primary_sig_comparison %>%
  
  dplyr::mutate(
    
    same_direction =
      sign(logFC_primary) ==
      sign(logFC_overlap)
    
  )


table(
  primary_sig_comparison$same_direction
)


cat(
  "Primary-significant genes tested in sensitivity:",
  nrow(primary_sig_comparison),
  "\n"
)


cat(
  "Same-direction genes:",
  sum(primary_sig_comparison$same_direction),
  "\n"
)


cat(
  "Directional concordance:",
  round(
    100 *
      mean(primary_sig_comparison$same_direction),
    1
  ),
  "%\n"
)
# ------------------------------------------------------------
# STEP 59: Effect-size concordance among the 95 comparable
#          primary-significant genes
# ------------------------------------------------------------

pearson_sig <- cor(
  primary_sig_comparison$logFC_primary,
  primary_sig_comparison$logFC_overlap,
  method = "pearson"
)


spearman_sig <- cor(
  primary_sig_comparison$logFC_primary,
  primary_sig_comparison$logFC_overlap,
  method = "spearman"
)


cat(
  "Pearson correlation:",
  round(pearson_sig, 4),
  "\n"
)


cat(
  "Spearman correlation:",
  round(spearman_sig, 4),
  "\n"
)
primary_sig_comparison <- primary_sig_comparison %>%
  
  dplyr::mutate(
    
    logFC_difference =
      logFC_overlap -
      logFC_primary,
    
    absolute_logFC_difference =
      abs(logFC_difference)
    
  )


summary(
  primary_sig_comparison$absolute_logFC_difference
)
# ------------------------------------------------------------
# STEP 61: Save overlapping-study sensitivity analysis
# ------------------------------------------------------------

readr::write_csv(
  de_overlap,
  "results/tables/VCM_age_DE_overlap15_study_adjusted.csv"
)


readr::write_csv(
  de_overlap_summary,
  "results/tables/VCM_age_DE_overlap15_summary.csv"
)


readr::write_csv(
  primary_sig_comparison,
  "results/tables/VCM_age_DE_primary_hits_sensitivity_comparison.csv"
)


readr::write_csv(
  primary_sig_not_tested,
  "results/tables/VCM_age_DE_primary_hits_not_tested_sensitivity.csv"
)


saveRDS(
  fit_overlap,
  "data/VCM_age_edgeR_QL_fit_overlap15.rds",
  compress = FALSE
)


saveRDS(
  qlf_overlap,
  "data/VCM_age_edgeR_QL_test_overlap15.rds",
  compress = FALSE
)
# ------------------------------------------------------------
# STEP 62: Create sensitivity-analysis summary
# ------------------------------------------------------------

robustness_summary <- tibble::tibble(
  
  primary_significant_genes = 111,
  
  primary_hits_tested_in_overlap =
    nrow(primary_sig_comparison),
  
  primary_hits_not_tested_in_overlap =
    nrow(primary_sig_not_tested),
  
  same_direction =
    sum(primary_sig_comparison$same_direction),
  
  directional_concordance_percent =
    100 * mean(primary_sig_comparison$same_direction),
  
  pearson_logFC =
    pearson_sig,
  
  spearman_logFC =
    spearman_sig,
  
  median_absolute_logFC_difference =
    median(
      primary_sig_comparison$absolute_logFC_difference
    ),
  
  maximum_absolute_logFC_difference =
    max(
      primary_sig_comparison$absolute_logFC_difference
    )
)


robustness_summary
readr::write_csv(
  robustness_summary,
  "results/tables/VCM_age_DE_sensitivity_robustness_summary.csv"
)
# ------------------------------------------------------------
# STEP 63: Load predefined biological panel
# ------------------------------------------------------------
#
# The biological panel was defined BEFORE inspecting the
# differential-expression results.
#
# This is important because it avoids selecting genes only
# because they happened to be significant.
#
# The panel includes:
#
#   - CYP enzymes
#   - EET / epoxide pathway
#   - HETE / omega-hydroxylation pathway
#   - arachidonic acid / COX / LOX genes
#   - estrogen signaling
#   - estrogen metabolism
#   - PPAR / PGC-1 metabolism
#   - AHR / xenobiotic signaling
#   - oxidative stress / redox genes
#
# If you already have the 194-gene vector saved in the
# current script/session, use it directly.
#
# Otherwise recreate/load it here before continuing.
ls()
# ============================================================
# TARGETED BIOLOGICAL PATHWAY ANALYSIS
# ============================================================
#
# PURPOSE:
#
# Interrogate a predefined panel of genes related to:
#
#   - CYP-mediated metabolism
#   - EET / epoxide metabolism
#   - HETE / omega-hydroxylation
#   - arachidonic acid / eicosanoid metabolism
#   - estrogen signaling
#   - steroid metabolism
#   - PPAR / PGC-1 signaling
#   - AHR / xenobiotic signaling
#   - oxidative stress / redox biology
#
# IMPORTANT:
#
# This panel is defined biologically, NOT by selecting genes
# because they were significant in the DE analysis.
#
# We compare:
#
#   PRIMARY MODEL:
#     all 24 donors
#     expression ~ DataSource + age_group
#
#   SENSITIVITY MODEL:
#     15 donors from studies containing both age groups
#
# Positive logFC:
#     higher expression in older females
#
# Negative logFC:
#     higher expression in younger females
#
# ============================================================


# ------------------------------------------------------------
# STEP 63: Identify all CYP genes represented in dataset
# ------------------------------------------------------------

all_dataset_genes <- rownames(pb_counts)


all_cyp_genes <- grep(
  "^CYP",
  all_dataset_genes,
  value = TRUE
)


cat(
  "Number of CYP genes represented in dataset:",
  length(all_cyp_genes),
  "\n"
)


# ------------------------------------------------------------
# STEP 64: Define additional pathway genes
# ------------------------------------------------------------

pathway_genes <- c(
  
  # ==========================================================
  # EET / epoxide metabolism
  # ==========================================================
  
  "EPHX2",
  "EPHX1",
  "POR",
  "CYB5A",
  "CYB5B",
  
  
  # ==========================================================
  # Arachidonic acid release / metabolism
  # ==========================================================
  
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
  
  "PTGS1",
  "PTGS2",
  
  "PTGES",
  "PTGES2",
  "PTGES3",
  
  "PTGDS",
  "HPGDS",
  "PTGIS",
  
  "TBXAS1",
  
  "HPGD",
  "PTGR1",
  "PTGR2",
  
  
  # ==========================================================
  # Lipoxygenase / leukotriene pathway
  # ==========================================================
  
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
  "GGT5",
  
  
  # ==========================================================
  # Eicosanoid receptors
  # ==========================================================
  
  "PTGER1",
  "PTGER2",
  "PTGER3",
  "PTGER4",
  
  "PTGDR",
  "PTGDR2",
  
  "PTGIR",
  
  "TBXA2R",
  
  "LTB4R",
  "LTB4R2",
  
  "CYSLTR1",
  "CYSLTR2",
  
  "OXER1",
  
  "GPR75",
  
  
  # ==========================================================
  # Estrogen signaling
  # ==========================================================
  
  "ESR1",
  "ESR2",
  "GPER1",
  
  "ESRRA",
  "ESRRB",
  "ESRRG",
  
  "NCOA1",
  "NCOA2",
  "NCOA3",
  
  "NCOR1",
  "NCOR2",
  
  
  # ==========================================================
  # Estrogen / steroid metabolism
  # ==========================================================
  
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
  
  "STS",
  "SULT1E1",
  
  "COMT",
  
  
  # ==========================================================
  # PPAR / PGC-1 metabolic signaling
  # ==========================================================
  
  "PPARA",
  "PPARD",
  "PPARG",
  
  "RXRA",
  "RXRB",
  "RXRG",
  
  "PPARGC1A",
  "PPARGC1B",
  
  "CD36",
  
  "FABP3",
  "FABP4",
  
  "CPT1A",
  "CPT1B",
  "CPT1C",
  
  "CPT2",
  
  "ACOX1",
  "ACOX2",
  
  "ACADL",
  "ACADM",
  "ACADVL",
  
  
  # ==========================================================
  # AHR signaling
  # ==========================================================
  
  "AHR",
  "ARNT",
  "ARNT2",
  
  "AHRR",
  
  "TIPARP",
  
  "NQO1",
  
  "ALDH3A1",
  
  
  # ==========================================================
  # Xenobiotic / nuclear receptor signaling
  # ==========================================================
  
  "NR1I2",
  "NR1I3",
  
  "NR1H4",
  "NR1H3",
  "NR1H2",
  
  "HNF4A",
  
  "RXRA",
  "RXRB",
  "RXRG",
  
  
  # ==========================================================
  # Oxidative stress / redox
  # ==========================================================
  
  "NFE2L2",
  
  "KEAP1",
  
  "NQO1",
  
  "HMOX1",
  "HMOX2",
  
  "SOD1",
  "SOD2",
  "SOD3",
  
  "CAT",
  
  "GPX1",
  "GPX2",
  "GPX3",
  "GPX4",
  
  "GCLC",
  "GCLM",
  
  "GSR",
  
  "TXN",
  
  "TXNRD1",
  "TXNRD2"
)


# ------------------------------------------------------------
# STEP 65: Combine all CYP genes with pathway genes
# ------------------------------------------------------------

biological_panel <- unique(
  c(
    all_cyp_genes,
    pathway_genes
  )
)


cat(
  "Total genes in reconstructed biological panel:",
  length(biological_panel),
  "\n"
)


# ------------------------------------------------------------
# STEP 66: Determine gene availability and testability
# ------------------------------------------------------------
#
# present_in_dataset:
#   gene exists in the processed feature matrix
#
# tested_primary:
#   gene passed filterByExpr in the 24-donor analysis
#
# tested_overlap:
#   gene passed filterByExpr in the 15-donor sensitivity
#   analysis

panel_availability <- tibble::tibble(
  
  gene =
    biological_panel,
  
  present_in_dataset =
    biological_panel %in% rownames(pb_counts),
  
  tested_primary =
    biological_panel %in% de_age$gene,
  
  tested_overlap =
    biological_panel %in% de_overlap$gene
)


# ------------------------------------------------------------
# STEP 67: Attach primary and sensitivity DE statistics
# ------------------------------------------------------------

target_results <- panel_availability %>%
  
  dplyr::left_join(
    
    de_age %>%
      dplyr::select(
        gene,
        logFC_primary = logFC,
        logCPM_primary = logCPM,
        PValue_primary = PValue,
        FDR_primary = FDR
      ),
    
    by = "gene"
    
  ) %>%
  
  dplyr::left_join(
    
    de_overlap %>%
      dplyr::select(
        gene,
        logFC_overlap = logFC,
        logCPM_overlap = logCPM,
        PValue_overlap = PValue,
        FDR_overlap = FDR
      ),
    
    by = "gene"
    
  )


# ------------------------------------------------------------
# STEP 68: Calculate directional consistency
# ------------------------------------------------------------
#
# same_direction = TRUE means the estimated older-vs-younger
# effect points in the same direction in both analyses.

target_results <- target_results %>%
  
  dplyr::mutate(
    
    same_direction = dplyr::case_when(
      
      tested_primary &
        tested_overlap ~
        
        sign(logFC_primary) ==
        sign(logFC_overlap),
      
      TRUE ~ NA
    )
    
  )


# ------------------------------------------------------------
# STEP 69: Classify level of evidence
# ------------------------------------------------------------

target_results <- target_results %>%
  
  dplyr::mutate(
    
    evidence_class = dplyr::case_when(
      
      tested_primary &
        FDR_primary < 0.05 &
        tested_overlap &
        same_direction ~
        
        "Primary FDR significant + directionally robust",
      
      tested_primary &
        FDR_primary < 0.05 &
        !tested_overlap ~
        
        "Primary FDR significant; not tested in sensitivity",
      
      tested_primary &
        PValue_primary < 0.05 &
        tested_overlap &
        same_direction ~
        
        "Nominal primary signal + directionally robust",
      
      tested_primary &
        PValue_primary < 0.05 ~
        
        "Nominal primary signal",
      
      tested_primary ~
        
        "Tested; no strong evidence",
      
      present_in_dataset &
        !tested_primary ~
        
        "Present but filtered from primary DE",
      
      TRUE ~
        
        "Absent from processed dataset"
    )
  )


# ------------------------------------------------------------
# STEP 70: Rank complete biological panel
# ------------------------------------------------------------

target_results_ranked <- target_results %>%
  
  dplyr::arrange(
    is.na(PValue_primary),
    PValue_primary
  )


# ------------------------------------------------------------
# STEP 71: Extract targeted genes significant at FDR < 0.05
# ------------------------------------------------------------

target_primary_sig <- target_results %>%
  
  dplyr::filter(
    tested_primary,
    FDR_primary < 0.05
  ) %>%
  
  dplyr::arrange(
    FDR_primary
  )


# ------------------------------------------------------------
# STEP 72: Extract nominal targeted signals
# ------------------------------------------------------------

target_nominal <- target_results %>%
  
  dplyr::filter(
    tested_primary,
    PValue_primary < 0.05
  ) %>%
  
  dplyr::arrange(
    PValue_primary
  )


# ------------------------------------------------------------
# STEP 73: Define core mechanistic genes
# ------------------------------------------------------------

core_mechanistic_genes <- c(
  
  # CYP / EET / sEH
  "CYP2J2",
  "EPHX2",
  "EPHX1",
  "POR",
  "CYB5A",
  
  # CYP / HETE
  "CYP4A11",
  "CYP4F11",
  "CYP4F12",
  "CYP4F22",
  "CYP2U1",
  
  # Estrogen / ERR
  "ESR1",
  "ESR2",
  "GPER1",
  "ESRRA",
  "ESRRB",
  "ESRRG",
  
  # PPAR / PGC1
  "PPARA",
  "PPARD",
  "PPARG",
  "PPARGC1A",
  "PPARGC1B",
  
  # AHR
  "AHR",
  "ARNT",
  "AHRR",
  "TIPARP",
  
  # Redox
  "NFE2L2",
  "KEAP1",
  "NQO1",
  "HMOX1",
  "SOD2"
)


core_results <- target_results %>%
  
  dplyr::filter(
    gene %in% core_mechanistic_genes
  ) %>%
  
  dplyr::mutate(
    
    gene = factor(
      gene,
      levels = core_mechanistic_genes
    )
    
  ) %>%
  
  dplyr::arrange(
    gene
  )


# ------------------------------------------------------------
# STEP 74: Print key diagnostics
# ------------------------------------------------------------

cat(
  "\n==============================\n"
)

cat(
  "TARGETED PANEL SUMMARY\n"
)

cat(
  "==============================\n"
)


cat(
  "CYP genes represented:",
  length(all_cyp_genes),
  "\n"
)


cat(
  "Total biological panel:",
  length(biological_panel),
  "\n"
)


cat(
  "Present in processed dataset:",
  sum(panel_availability$present_in_dataset),
  "\n"
)


cat(
  "Tested in primary analysis:",
  sum(panel_availability$tested_primary),
  "\n"
)


cat(
  "Tested in sensitivity analysis:",
  sum(panel_availability$tested_overlap),
  "\n"
)


cat(
  "Targeted genes with primary P < 0.05:",
  nrow(target_nominal),
  "\n"
)


cat(
  "Targeted genes with primary FDR < 0.05:",
  nrow(target_primary_sig),
  "\n"
)


# ------------------------------------------------------------
# STEP 75: Print FDR-significant targeted genes
# ------------------------------------------------------------

cat(
  "\n==============================\n"
)

cat(
  "TARGETED FDR-SIGNIFICANT GENES\n"
)

cat(
  "==============================\n"
)


target_primary_sig %>%
  
  dplyr::select(
    gene,
    logFC_primary,
    FDR_primary,
    logFC_overlap,
    PValue_overlap,
    FDR_overlap,
    same_direction,
    evidence_class
  ) %>%
  
  tibble::as_tibble() %>%
  
  print(
    n = Inf
  )


# ------------------------------------------------------------
# STEP 76: Print core mechanistic genes
# ------------------------------------------------------------

cat(
  "\n==============================\n"
)

cat(
  "CORE MECHANISTIC GENES\n"
)

cat(
  "==============================\n"
)


core_results %>%
  
  dplyr::select(
    gene,
    present_in_dataset,
    tested_primary,
    logFC_primary,
    PValue_primary,
    FDR_primary,
    tested_overlap,
    logFC_overlap,
    PValue_overlap,
    FDR_overlap,
    same_direction,
    evidence_class
  ) %>%
  
  tibble::as_tibble() %>%
  
  print(
    n = Inf
  )


# ------------------------------------------------------------
# STEP 77: Save targeted pathway tables
# ------------------------------------------------------------

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)


readr::write_csv(
  target_results_ranked,
  "results/tables/VCM_age_targeted_biological_panel_all.csv"
)


readr::write_csv(
  target_nominal,
  "results/tables/VCM_age_targeted_biological_panel_nominal.csv"
)


readr::write_csv(
  target_primary_sig,
  "results/tables/VCM_age_targeted_biological_panel_FDRsig.csv"
)


readr::write_csv(
  core_results,
  "results/tables/VCM_age_core_mechanistic_genes.csv"
)


# ------------------------------------------------------------
# STEP 78: Save biological panel itself
# ------------------------------------------------------------

saveRDS(
  biological_panel,
  "data/VCM_predefined_biological_panel.rds"
)


# ------------------------------------------------------------
# FINAL MESSAGE
# ------------------------------------------------------------

cat(
  "\nTargeted biological pathway analysis complete.\n"
)
# ============================================================
# TARGETED PATHWAY RESULTS:
# FIGURES + DETAILED OUTPUT
# ============================================================
#
# This section:
#
#   1. Prints all 21 nominal targeted genes
#   2. Examines sensitivity-analysis consistency
#   3. Creates a focused core-gene effect plot
#   4. Creates a targeted nominal-gene effect plot
#   5. Saves publication/presentation-ready figures
#
# IMPORTANT:
#
# Positive log2FC = higher expression in OLDER females
# Negative log2FC = higher expression in YOUNGER females
#
# ============================================================


# ------------------------------------------------------------
# STEP 79: Print all targeted genes with primary P < 0.05
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "TARGETED GENES WITH PRIMARY P < 0.05\n",
  "============================================\n"
)


target_nominal %>%
  
  dplyr::select(
    gene,
    logFC_primary,
    logCPM_primary,
    PValue_primary,
    FDR_primary,
    tested_overlap,
    logFC_overlap,
    PValue_overlap,
    FDR_overlap,
    same_direction,
    evidence_class
  ) %>%
  
  tibble::as_tibble() %>%
  
  print(
    n = Inf,
    width = Inf
  )


# ------------------------------------------------------------
# STEP 80: Quantify robustness of the targeted nominal genes
# ------------------------------------------------------------

target_nominal_testable <- target_nominal %>%
  
  dplyr::filter(
    tested_overlap
  )


cat(
  "\n============================================\n",
  "TARGETED NOMINAL-GENE ROBUSTNESS\n",
  "============================================\n"
)


cat(
  "Primary nominal targeted genes:",
  nrow(target_nominal),
  "\n"
)


cat(
  "Also tested in sensitivity:",
  nrow(target_nominal_testable),
  "\n"
)


cat(
  "Same direction:",
  sum(
    target_nominal_testable$same_direction,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Directional concordance:",
  round(
    100 *
      mean(
        target_nominal_testable$same_direction,
        na.rm = TRUE
      ),
    1
  ),
  "%\n"
)


# ------------------------------------------------------------
# STEP 81: Correlation for targeted nominal genes
# ------------------------------------------------------------

if (nrow(target_nominal_testable) >= 3) {
  
  target_nominal_pearson <- cor(
    target_nominal_testable$logFC_primary,
    target_nominal_testable$logFC_overlap,
    method = "pearson"
  )
  
  target_nominal_spearman <- cor(
    target_nominal_testable$logFC_primary,
    target_nominal_testable$logFC_overlap,
    method = "spearman"
  )
  
  
  cat(
    "Pearson logFC correlation:",
    round(target_nominal_pearson, 4),
    "\n"
  )
  
  
  cat(
    "Spearman logFC correlation:",
    round(target_nominal_spearman, 4),
    "\n"
  )
  
}


# ------------------------------------------------------------
# STEP 82: Create plotting table for core mechanistic genes
# ------------------------------------------------------------

core_plot_data <- core_results %>%
  
  dplyr::filter(
    tested_primary
  ) %>%
  
  dplyr::mutate(
    
    gene =
      factor(
        as.character(gene),
        levels = rev(
          core_mechanistic_genes[
            core_mechanistic_genes %in%
              as.character(gene)
          ]
        )
      ),
    
    significance =
      dplyr::case_when(
        
        FDR_primary < 0.05 ~
          "FDR < 0.05",
        
        PValue_primary < 0.05 ~
          "Nominal P < 0.05",
        
        TRUE ~
          "Not significant"
      )
  )


# ------------------------------------------------------------
# STEP 83: Core mechanistic-gene forest-style plot
# ------------------------------------------------------------

p_core_effects <- ggplot(
  core_plot_data,
  aes(
    x = logFC_primary,
    y = gene
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(
      shape = significance
    ),
    size = 3
  ) +
  
  labs(
    
    title =
      "Age-associated expression of core metabolic genes",
    
    subtitle =
      "Female LV ventricular cardiomyocytes; study-adjusted pseudobulk model",
    
    x =
      "log2 fold-change (older vs younger)",
    
    y =
      NULL,
    
    shape =
      "Primary analysis"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 15
      ),
    
    axis.text.y =
      element_text(
        face = "italic"
      ),
    
    legend.position =
      "bottom"
  )


print(
  p_core_effects
)


# ------------------------------------------------------------
# STEP 84: Prepare the 21 nominal targeted genes for plotting
# ------------------------------------------------------------

target_nominal_plot <- target_nominal %>%
  
  dplyr::mutate(
    
    significance =
      dplyr::if_else(
        FDR_primary < 0.05,
        "FDR < 0.05",
        "Nominal P < 0.05"
      ),
    
    gene =
      factor(
        gene,
        levels =
          rev(
            gene[
              order(logFC_primary)
            ]
          )
      )
  )


# ------------------------------------------------------------
# STEP 85: Plot all nominal targeted signals
# ------------------------------------------------------------

p_target_nominal <- ggplot(
  target_nominal_plot,
  aes(
    x = logFC_primary,
    y = gene
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_point(
    aes(
      shape = significance
    ),
    size = 3.2
  ) +
  
  labs(
    
    title =
      "Targeted pathway genes associated with age",
    
    subtitle =
      "Genes with primary P < 0.05; positive values indicate higher expression in older females",
    
    x =
      "log2 fold-change (older vs younger)",
    
    y =
      NULL,
    
    shape =
      "Evidence"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 15
      ),
    
    axis.text.y =
      element_text(
        face = "italic"
      ),
    
    legend.position =
      "bottom"
  )


print(
  p_target_nominal
)


# ------------------------------------------------------------
# STEP 86: Primary vs sensitivity plot for targeted nominal
#          genes
# ------------------------------------------------------------

p_target_robustness <- ggplot(
  
  target_nominal_testable,
  
  aes(
    x = logFC_primary,
    y = logFC_overlap
  )
  
) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  
  geom_point(
    size = 3
  ) +
  
  ggrepel::geom_text_repel(
    
    aes(
      label = gene
    ),
    
    size = 3.5,
    
    max.overlaps = Inf
  ) +
  
  labs(
    
    title =
      "Robustness of targeted age-associated signals",
    
    subtitle =
      "Primary versus overlapping-study sensitivity analysis",
    
    x =
      "log2FC: all 24 donors",
    
    y =
      "log2FC: overlapping studies (15 donors)"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 15
      )
  )


print(
  p_target_robustness
)


# ------------------------------------------------------------
# STEP 87: Save figures
# ------------------------------------------------------------

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)


ggplot2::ggsave(
  
  filename =
    "results/figures/VCM_age_core_mechanistic_genes.pdf",
  
  plot =
    p_core_effects,
  
  width = 8,
  height = 8
)


ggplot2::ggsave(
  
  filename =
    "results/figures/VCM_age_core_mechanistic_genes.png",
  
  plot =
    p_core_effects,
  
  width = 8,
  height = 8,
  
  dpi = 300
)


ggplot2::ggsave(
  
  filename =
    "results/figures/VCM_age_targeted_nominal_genes.pdf",
  
  plot =
    p_target_nominal,
  
  width = 8,
  height = 7
)


ggplot2::ggsave(
  
  filename =
    "results/figures/VCM_age_targeted_nominal_genes.png",
  
  plot =
    p_target_nominal,
  
  width = 8,
  height = 7,
  
  dpi = 300
)


ggplot2::ggsave(
  
  filename =
    "results/figures/VCM_age_targeted_sensitivity_concordance.pdf",
  
  plot =
    p_target_robustness,
  
  width = 8,
  height = 7
)


ggplot2::ggsave(
  
  filename =
    "results/figures/VCM_age_targeted_sensitivity_concordance.png",
  
  plot =
    p_target_robustness,
  
  width = 8,
  height = 7,
  
  dpi = 300
)


# ------------------------------------------------------------
# STEP 88: Save targeted robustness summary
# ------------------------------------------------------------

target_robustness_summary <- tibble::tibble(
  
  targeted_nominal_primary =
    nrow(target_nominal),
  
  targeted_nominal_tested_overlap =
    nrow(target_nominal_testable),
  
  same_direction =
    sum(
      target_nominal_testable$same_direction,
      na.rm = TRUE
    ),
  
  directional_concordance_percent =
    100 *
    mean(
      target_nominal_testable$same_direction,
      na.rm = TRUE
    ),
  
  pearson_logFC =
    target_nominal_pearson,
  
  spearman_logFC =
    target_nominal_spearman
)


print(
  target_robustness_summary
)


readr::write_csv(
  
  target_robustness_summary,
  
  "results/tables/VCM_age_targeted_robustness_summary.csv"
)


# ------------------------------------------------------------
# FINAL MESSAGE
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "TARGETED FIGURE ANALYSIS COMPLETE\n",
  "============================================\n"
)

# ============================================================
# STEP 89: CHECK ENRICHMENT PACKAGES
# ============================================================
#
# We will perform unbiased genome-wide pathway enrichment
# using the complete ranked edgeR result.
#
# First determine which enrichment packages are available
# in the current R installation.
#
# No packages are installed or modified by this block.
# ============================================================

enrichment_packages <- c(
  "fgsea",
  "msigdbr",
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA"
)


package_status <- tibble::tibble(
  
  package = enrichment_packages,
  
  installed = vapply(
    enrichment_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
  
)


print(
  package_status,
  n = Inf
)


# ------------------------------------------------------------
# Also confirm the edgeR result contains the statistic needed
# for ranking.
# ------------------------------------------------------------

cat(
  "\nColumns in de_age:\n"
)

print(
  colnames(de_age)
)


cat(
  "\nNumber of genes in primary DE result:",
  nrow(de_age),
  "\n"
)


cat(
  "\nF-statistic summary:\n"
)

print(
  summary(de_age$F)
)
# ============================================================
# UNBIASED GENOME-WIDE PATHWAY ENRICHMENT
# ============================================================
#
# PURPOSE:
#
# Use ALL genes tested in the primary edgeR pseudobulk model,
# rather than selecting only FDR-significant genes.
#
# We will run preranked GSEA against:
#
#   1. MSigDB Hallmark pathways
#   2. Reactome pathways
#   3. Gene Ontology Biological Process (GO BP)
#
# PRIMARY DE MODEL:
#
#   all 24 female LV donors
#   expression ~ DataSource + age_group
#
# AGE COEFFICIENT:
#
#   Older_55plus versus Younger_20_45
#
# Therefore:
#
#   Positive ranking statistic / positive NES
#       = enriched in OLDER female VCMs
#
#   Negative ranking statistic / negative NES
#       = enriched in YOUNGER female VCMs
#
# ============================================================



# ============================================================
# STEP 89A: INSTALL REQUIRED PACKAGES IF MISSING
# ============================================================

# ------------------------------------------------------------
# Install BiocManager if necessary
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  
  install.packages(
    "BiocManager"
  )
  
}


# ------------------------------------------------------------
# Install fgsea from Bioconductor if necessary
# ------------------------------------------------------------

if (!requireNamespace("fgsea", quietly = TRUE)) {
  
  BiocManager::install(
    "fgsea",
    ask = FALSE,
    update = FALSE
  )
  
}


# ------------------------------------------------------------
# Install msigdbr from CRAN if necessary
# ------------------------------------------------------------

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  
  install.packages(
    "msigdbr"
  )
  
}


# ------------------------------------------------------------
# Confirm both packages are now available
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "PACKAGE CHECK\n",
  "============================================\n"
)


cat(
  "fgsea installed:",
  requireNamespace(
    "fgsea",
    quietly = TRUE
  ),
  "\n"
)


cat(
  "msigdbr installed:",
  requireNamespace(
    "msigdbr",
    quietly = TRUE
  ),
  "\n"
)


# Stop immediately if installation did not succeed.

stopifnot(
  requireNamespace("fgsea", quietly = TRUE),
  requireNamespace("msigdbr", quietly = TRUE)
)



# ============================================================
# STEP 90: BUILD THE GENOME-WIDE RANKING STATISTIC
# ============================================================
#
# edgeR's quasi-likelihood F statistic is always positive.
#
# To use it for directional enrichment, multiply sqrt(F)
# by the sign of the age-associated log2FC.
#
# This creates a t-like signed statistic:
#
#       sign(logFC) * sqrt(F)
#
# Positive:
#       higher in OLDER donors
#
# Negative:
#       higher in YOUNGER donors
#
# We use all genes retained by filterByExpr.
# ============================================================


gsea_rank_table <- de_age %>%
  
  dplyr::filter(
    !is.na(logFC),
    !is.na(F),
    is.finite(logFC),
    is.finite(F)
  ) %>%
  
  dplyr::mutate(
    
    rank_stat =
      sign(logFC) *
      sqrt(F)
    
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(rank_stat)
  )


cat(
  "\n============================================\n",
  "GSEA RANKING TABLE\n",
  "============================================\n"
)


cat(
  "Genes available for ranking:",
  nrow(gsea_rank_table),
  "\n"
)


print(
  summary(
    gsea_rank_table$rank_stat
  )
)


# ------------------------------------------------------------
# Create named numeric vector required by fgsea
# ------------------------------------------------------------

gsea_ranks <- gsea_rank_table$rank_stat


names(gsea_ranks) <-
  gsea_rank_table$gene


# Sort from most older-associated to most younger-associated

gsea_ranks <- sort(
  gsea_ranks,
  decreasing = TRUE
)


# ------------------------------------------------------------
# Verify ranking
# ------------------------------------------------------------

cat(
  "\nTop 10 genes toward OLDER females:\n"
)


print(
  head(
    gsea_ranks,
    10
  )
)


cat(
  "\nTop 10 genes toward YOUNGER females:\n"
)


print(
  tail(
    gsea_ranks,
    10
  )
)



# ============================================================
# STEP 91: DOWNLOAD/LOAD HALLMARK GENE SETS FROM MSigDB
# ============================================================


hallmark_msig <- msigdbr::msigdbr(
  
  species = "Homo sapiens",
  
  collection = "H"
  
)


cat(
  "\n============================================\n",
  "HALLMARK COLLECTION\n",
  "============================================\n"
)


cat(
  "Hallmark gene-set rows:",
  nrow(hallmark_msig),
  "\n"
)


cat(
  "Hallmark pathways:",
  dplyr::n_distinct(
    hallmark_msig$gs_name
  ),
  "\n"
)


# Convert MSigDB table into fgsea pathway list

hallmark_pathways <- split(
  
  hallmark_msig$gene_symbol,
  
  hallmark_msig$gs_name
  
)



# ============================================================
# STEP 92: RUN HALLMARK GSEA
# ============================================================


set.seed(
  12345
)


fgsea_hallmark <- fgsea::fgseaMultilevel(
  
  pathways =
    hallmark_pathways,
  
  stats =
    gsea_ranks,
  
  minSize =
    10,
  
  maxSize =
    500,
  
  eps =
    0
  
)


fgsea_hallmark <- fgsea_hallmark %>%
  
  tibble::as_tibble() %>%
  
  dplyr::arrange(
    padj
  )


cat(
  "\n============================================\n",
  "HALLMARK GSEA SUMMARY\n",
  "============================================\n"
)


cat(
  "Pathways tested:",
  nrow(fgsea_hallmark),
  "\n"
)


cat(
  "FDR < 0.05:",
  sum(
    fgsea_hallmark$padj < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "FDR < 0.05 enriched in older:",
  sum(
    fgsea_hallmark$padj < 0.05 &
      fgsea_hallmark$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "FDR < 0.05 enriched in younger:",
  sum(
    fgsea_hallmark$padj < 0.05 &
      fgsea_hallmark$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nTop Hallmark pathways:\n"
)


fgsea_hallmark %>%
  
  dplyr::select(
    pathway,
    size,
    NES,
    pval,
    padj
  ) %>%
  
  dplyr::slice_head(
    n = 25
  ) %>%
  
  print(
    n = 25,
    width = Inf
  )



# ============================================================
# STEP 93: LOAD REACTOME GENE SETS
# ============================================================


reactome_msig <- msigdbr::msigdbr(
  
  species = "Homo sapiens",
  
  collection = "C2",
  
  subcollection = "CP:REACTOME"
  
)


cat(
  "\n============================================\n",
  "REACTOME COLLECTION\n",
  "============================================\n"
)


cat(
  "Reactome pathways:",
  dplyr::n_distinct(
    reactome_msig$gs_name
  ),
  "\n"
)


reactome_pathways <- split(
  
  reactome_msig$gene_symbol,
  
  reactome_msig$gs_name
  
)



# ============================================================
# STEP 94: RUN REACTOME GSEA
# ============================================================


set.seed(
  12345
)


fgsea_reactome <- fgsea::fgseaMultilevel(
  
  pathways =
    reactome_pathways,
  
  stats =
    gsea_ranks,
  
  minSize =
    10,
  
  maxSize =
    500,
  
  eps =
    0
  
)


fgsea_reactome <- fgsea_reactome %>%
  
  tibble::as_tibble() %>%
  
  dplyr::arrange(
    padj
  )


cat(
  "\n============================================\n",
  "REACTOME GSEA SUMMARY\n",
  "============================================\n"
)


cat(
  "Pathways tested:",
  nrow(fgsea_reactome),
  "\n"
)


cat(
  "FDR < 0.05:",
  sum(
    fgsea_reactome$padj < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "FDR < 0.05 enriched in older:",
  sum(
    fgsea_reactome$padj < 0.05 &
      fgsea_reactome$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "FDR < 0.05 enriched in younger:",
  sum(
    fgsea_reactome$padj < 0.05 &
      fgsea_reactome$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nTop Reactome pathways:\n"
)


fgsea_reactome %>%
  
  dplyr::select(
    pathway,
    size,
    NES,
    pval,
    padj
  ) %>%
  
  dplyr::slice_head(
    n = 30
  ) %>%
  
  print(
    n = 30,
    width = Inf
  )



# ============================================================
# STEP 95: LOAD GO BIOLOGICAL PROCESS GENE SETS
# ============================================================


gobp_msig <- msigdbr::msigdbr(
  
  species = "Homo sapiens",
  
  collection = "C5",
  
  subcollection = "GO:BP"
  
)


cat(
  "\n============================================\n",
  "GO BIOLOGICAL PROCESS COLLECTION\n",
  "============================================\n"
)


cat(
  "GO BP pathways:",
  dplyr::n_distinct(
    gobp_msig$gs_name
  ),
  "\n"
)


gobp_pathways <- split(
  
  gobp_msig$gene_symbol,
  
  gobp_msig$gs_name
  
)



# ============================================================
# STEP 96: RUN GO BIOLOGICAL PROCESS GSEA
# ============================================================


set.seed(
  12345
)


fgsea_gobp <- fgsea::fgseaMultilevel(
  
  pathways =
    gobp_pathways,
  
  stats =
    gsea_ranks,
  
  minSize =
    10,
  
  maxSize =
    500,
  
  eps =
    0
  
)


fgsea_gobp <- fgsea_gobp %>%
  
  tibble::as_tibble() %>%
  
  dplyr::arrange(
    padj
  )


cat(
  "\n============================================\n",
  "GO BP GSEA SUMMARY\n",
  "============================================\n"
)


cat(
  "Pathways tested:",
  nrow(fgsea_gobp),
  "\n"
)


cat(
  "FDR < 0.05:",
  sum(
    fgsea_gobp$padj < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "FDR < 0.05 enriched in older:",
  sum(
    fgsea_gobp$padj < 0.05 &
      fgsea_gobp$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "FDR < 0.05 enriched in younger:",
  sum(
    fgsea_gobp$padj < 0.05 &
      fgsea_gobp$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nTop GO BP pathways:\n"
)


fgsea_gobp %>%
  
  dplyr::select(
    pathway,
    size,
    NES,
    pval,
    padj
  ) %>%
  
  dplyr::slice_head(
    n = 40
  ) %>%
  
  print(
    n = 40,
    width = Inf
  )



# ============================================================
# STEP 97: ADD ENRICHMENT DIRECTION LABELS
# ============================================================


fgsea_hallmark <- fgsea_hallmark %>%
  
  dplyr::mutate(
    
    direction =
      dplyr::if_else(
        NES > 0,
        "Older",
        "Younger"
      )
    
  )


fgsea_reactome <- fgsea_reactome %>%
  
  dplyr::mutate(
    
    direction =
      dplyr::if_else(
        NES > 0,
        "Older",
        "Younger"
      )
    
  )


fgsea_gobp <- fgsea_gobp %>%
  
  dplyr::mutate(
    
    direction =
      dplyr::if_else(
        NES > 0,
        "Older",
        "Younger"
      )
    
  )



# ============================================================
# STEP 98: EXTRACT SIGNIFICANT PATHWAYS
# ============================================================


hallmark_sig <- fgsea_hallmark %>%
  
  dplyr::filter(
    padj < 0.05
  ) %>%
  
  dplyr::arrange(
    padj
  )


reactome_sig <- fgsea_reactome %>%
  
  dplyr::filter(
    padj < 0.05
  ) %>%
  
  dplyr::arrange(
    padj
  )


gobp_sig <- fgsea_gobp %>%
  
  dplyr::filter(
    padj < 0.05
  ) %>%
  
  dplyr::arrange(
    padj
  )



# ============================================================
# STEP 99: SEARCH ENRICHMENT RESULTS FOR OUR BIOLOGICAL THEMES
# ============================================================
#
# This does NOT determine significance.
#
# It simply lets us inspect whether pathways relevant to our
# original biological question appear anywhere in the
# unbiased genome-wide enrichment results.
#
# ============================================================


pathway_keywords <- paste(
  
  c(
    "CYTOCHROME",
    "ARACHID",
    "EICOS",
    "PROSTAGLAND",
    "LEUKOTRIENE",
    "LIPOXYGEN",
    "FATTY_ACID",
    "LIPID",
    "ESTROGEN",
    "PPAR",
    "PEROXIS",
    "XENOBIOTIC",
    "OXIDATIVE",
    "REACTIVE_OXYGEN",
    "GLUTATHIONE",
    "MITOCHON",
    "CARDIAC",
    "MUSCLE"
  ),
  
  collapse = "|"
  
)


hallmark_theme_hits <- fgsea_hallmark %>%
  
  dplyr::filter(
    grepl(
      pathway_keywords,
      pathway,
      ignore.case = TRUE
    )
  ) %>%
  
  dplyr::arrange(
    padj
  )


reactome_theme_hits <- fgsea_reactome %>%
  
  dplyr::filter(
    grepl(
      pathway_keywords,
      pathway,
      ignore.case = TRUE
    )
  ) %>%
  
  dplyr::arrange(
    padj
  )


gobp_theme_hits <- fgsea_gobp %>%
  
  dplyr::filter(
    grepl(
      pathway_keywords,
      pathway,
      ignore.case = TRUE
    )
  ) %>%
  
  dplyr::arrange(
    padj
  )


cat(
  "\n============================================\n",
  "BIOLOGICALLY RELEVANT HALLMARK RESULTS\n",
  "============================================\n"
)


hallmark_theme_hits %>%
  
  dplyr::select(
    pathway,
    NES,
    pval,
    padj,
    direction
  ) %>%
  
  print(
    n = Inf,
    width = Inf
  )


cat(
  "\n============================================\n",
  "BIOLOGICALLY RELEVANT REACTOME RESULTS\n",
  "============================================\n"
)


reactome_theme_hits %>%
  
  dplyr::select(
    pathway,
    NES,
    pval,
    padj,
    direction
  ) %>%
  
  print(
    n = Inf,
    width = Inf
  )


cat(
  "\n============================================\n",
  "BIOLOGICALLY RELEVANT GO BP RESULTS\n",
  "============================================\n"
)


gobp_theme_hits %>%
  
  dplyr::select(
    pathway,
    NES,
    pval,
    padj,
    direction
  ) %>%
  
  print(
    n = Inf,
    width = Inf
  )



# ============================================================
# STEP 100: MAKE A CLEAN HALLMARK ENRICHMENT FIGURE
# ============================================================
#
# Hallmark is usually the cleanest collection for a
# presentation because it contains only ~50 non-redundant
# biological programs.
#
# We plot the strongest pathways ranked by FDR.
# ============================================================


hallmark_plot_data <- fgsea_hallmark %>%
  
  dplyr::filter(
    !is.na(padj)
  ) %>%
  
  dplyr::arrange(
    padj
  ) %>%
  
  dplyr::slice_head(
    n = 20
  ) %>%
  
  dplyr::mutate(
    
    pathway_clean =
      gsub(
        "^HALLMARK_",
        "",
        pathway
      ),
    
    pathway_clean =
      gsub(
        "_",
        " ",
        pathway_clean
      ),
    
    pathway_clean =
      factor(
        pathway_clean,
        levels =
          rev(pathway_clean)
      )
    
  )


p_hallmark_gsea <- ggplot2::ggplot(
  
  hallmark_plot_data,
  
  ggplot2::aes(
    x = NES,
    y = pathway_clean
  )
  
) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  ggplot2::geom_point(
    
    ggplot2::aes(
      size = -log10(padj),
      shape = direction
    )
    
  ) +
  
  ggplot2::labs(
    
    title =
      "Genome-wide Hallmark pathway enrichment",
    
    subtitle =
      "Female LV ventricular cardiomyocytes: older versus younger donors",
    
    x =
      "Normalized enrichment score (NES)",
    
    y =
      NULL,
    
    size =
      "-log10(FDR)",
    
    shape =
      "Enriched toward"
    
  ) +
  
  ggplot2::theme_classic(
    base_size = 13
  ) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    
    legend.position =
      "bottom"
    
  )


print(
  p_hallmark_gsea
)



# ============================================================
# STEP 101: SAVE ALL ENRICHMENT RESULTS
# ============================================================


dir.create(
  "results/enrichment",
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Save ranking table
# ------------------------------------------------------------

readr::write_csv(
  
  gsea_rank_table,
  
  "results/enrichment/VCM_age_GSEA_ranked_genes.csv"
  
)


# ------------------------------------------------------------
# Save complete enrichment tables
# ------------------------------------------------------------

readr::write_csv(
  
  fgsea_hallmark,
  
  "results/enrichment/VCM_age_GSEA_Hallmark_all.csv"
  
)


readr::write_csv(
  
  fgsea_reactome,
  
  "results/enrichment/VCM_age_GSEA_Reactome_all.csv"
  
)


readr::write_csv(
  
  fgsea_gobp,
  
  "results/enrichment/VCM_age_GSEA_GO_BP_all.csv"
  
)


# ------------------------------------------------------------
# Save FDR-significant pathway tables
# ------------------------------------------------------------

readr::write_csv(
  
  hallmark_sig,
  
  "results/enrichment/VCM_age_GSEA_Hallmark_FDR05.csv"
  
)


readr::write_csv(
  
  reactome_sig,
  
  "results/enrichment/VCM_age_GSEA_Reactome_FDR05.csv"
  
)


readr::write_csv(
  
  gobp_sig,
  
  "results/enrichment/VCM_age_GSEA_GO_BP_FDR05.csv"
  
)


# ------------------------------------------------------------
# Save biologically relevant keyword results
# ------------------------------------------------------------

readr::write_csv(
  
  hallmark_theme_hits,
  
  "results/enrichment/VCM_age_GSEA_Hallmark_target_themes.csv"
  
)


readr::write_csv(
  
  reactome_theme_hits,
  
  "results/enrichment/VCM_age_GSEA_Reactome_target_themes.csv"
  
)


readr::write_csv(
  
  gobp_theme_hits,
  
  "results/enrichment/VCM_age_GSEA_GO_BP_target_themes.csv"
  
)


# ------------------------------------------------------------
# Save Hallmark figure
# ------------------------------------------------------------

ggplot2::ggsave(
  
  filename =
    "results/figures/VCM_age_GSEA_Hallmark_top20.pdf",
  
  plot =
    p_hallmark_gsea,
  
  width =
    10,
  
  height =
    8
  
)


ggplot2::ggsave(
  
  filename =
    "results/figures/VCM_age_GSEA_Hallmark_top20.png",
  
  plot =
    p_hallmark_gsea,
  
  width =
    10,
  
  height =
    8,
  
  dpi =
    300
  
)



# ============================================================
# STEP 102: FINAL COMPACT SUMMARY
# ============================================================


cat(
  "\n\n============================================\n",
  "GENOME-WIDE GSEA COMPLETE\n",
  "============================================\n"
)


cat(
  "\nHALLMARK\n"
)


cat(
  "FDR-significant pathways:",
  nrow(hallmark_sig),
  "\n"
)


cat(
  "Older enriched:",
  sum(
    hallmark_sig$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Younger enriched:",
  sum(
    hallmark_sig$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nREACTOME\n"
)


cat(
  "FDR-significant pathways:",
  nrow(reactome_sig),
  "\n"
)


cat(
  "Older enriched:",
  sum(
    reactome_sig$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Younger enriched:",
  sum(
    reactome_sig$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nGO BIOLOGICAL PROCESS\n"
)


cat(
  "FDR-significant pathways:",
  nrow(gobp_sig),
  "\n"
)


cat(
  "Older enriched:",
  sum(
    gobp_sig$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Younger enriched:",
  sum(
    gobp_sig$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\n============================================\n",
  "DONE\n",
  "============================================\n"
)

# ============================================================
# FIX GSEA RANKING: protect against tiny negative F values
# ============================================================

gsea_rank_table <- de_age %>%
  
  dplyr::filter(
    !is.na(logFC),
    !is.na(F),
    is.finite(logFC),
    is.finite(F)
  ) %>%
  
  dplyr::mutate(
    
    # edgeR F statistics should theoretically be >= 0.
    # pmax() prevents tiny numerical negative values from
    # generating NaN when square-rooted.
    
    F_nonnegative =
      pmax(F, 0),
    
    rank_stat =
      sign(logFC) *
      sqrt(F_nonnegative)
    
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(rank_stat)
  )


cat(
  "NA ranking statistics:",
  sum(is.na(gsea_rank_table$rank_stat)),
  "\n"
)

cat(
  "Negative raw F values:",
  sum(de_age$F < 0, na.rm = TRUE),
  "\n"
)


gsea_ranks <- gsea_rank_table$rank_stat
names(gsea_ranks) <- gsea_rank_table$gene

gsea_ranks <- sort(
  gsea_ranks,
  decreasing = TRUE
)
# ============================================================
# FINAL GSEA RERUN AFTER FIXING NEGATIVE F VALUES
# ============================================================
#
# RUN THIS BLOCK AFTER:
#
#   gsea_rank_table
#   gsea_ranks
#
# have been recreated using:
#
#   F_nonnegative = pmax(F, 0)
#
# This reruns:
#
#   1. Hallmark GSEA
#   2. Reactome GSEA
#   3. GO Biological Process GSEA
#   4. Direction labels
#   5. FDR-significant result tables
#   6. Target-theme extraction
#   7. Clean Hallmark figure
#   8. Final saved output files
#
# Positive NES = enriched toward OLDER females
# Negative NES = enriched toward YOUNGER females
#
# ============================================================


# ------------------------------------------------------------
# STEP 1: Sanity-check corrected ranking vector
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "CORRECTED GSEA RANK CHECK\n",
  "============================================\n"
)

cat(
  "Genes in ranking vector:",
  length(gsea_ranks),
  "\n"
)

cat(
  "NA values:",
  sum(is.na(gsea_ranks)),
  "\n"
)

cat(
  "Infinite values:",
  sum(!is.finite(gsea_ranks)),
  "\n"
)

cat(
  "Duplicated gene names:",
  sum(duplicated(names(gsea_ranks))),
  "\n"
)


# We want:
#
# NA values = 0
# Infinite values = 0


# ============================================================
# STEP 2: HALLMARK
# ============================================================

# Recreate Hallmark pathway collection if necessary

if (!exists("hallmark_pathways")) {
  
  hallmark_msig <- msigdbr::msigdbr(
    species = "Homo sapiens",
    collection = "H"
  )
  
  hallmark_pathways <- split(
    hallmark_msig$gene_symbol,
    hallmark_msig$gs_name
  )
}


# ------------------------------------------------------------
# Run Hallmark GSEA
# ------------------------------------------------------------

set.seed(12345)

fgsea_hallmark <- fgsea::fgseaMultilevel(
  
  pathways = hallmark_pathways,
  
  stats = gsea_ranks,
  
  minSize = 10,
  
  maxSize = 500,
  
  eps = 0
)


fgsea_hallmark <- fgsea_hallmark %>%
  
  tibble::as_tibble() %>%
  
  dplyr::mutate(
    
    direction =
      dplyr::if_else(
        NES > 0,
        "Older",
        "Younger"
      )
    
  ) %>%
  
  dplyr::arrange(
    padj
  )


# ------------------------------------------------------------
# Significant Hallmark pathways
# ------------------------------------------------------------

hallmark_sig <- fgsea_hallmark %>%
  
  dplyr::filter(
    padj < 0.05
  ) %>%
  
  dplyr::arrange(
    padj
  )


# ------------------------------------------------------------
# Print Hallmark summary
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "FINAL HALLMARK GSEA SUMMARY\n",
  "============================================\n"
)

cat(
  "Pathways tested:",
  nrow(fgsea_hallmark),
  "\n"
)

cat(
  "FDR < 0.05:",
  nrow(hallmark_sig),
  "\n"
)

cat(
  "Older enriched:",
  sum(
    hallmark_sig$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Younger enriched:",
  sum(
    hallmark_sig$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nTop 20 Hallmark pathways:\n"
)

fgsea_hallmark %>%
  
  dplyr::select(
    pathway,
    size,
    NES,
    pval,
    padj,
    direction
  ) %>%
  
  dplyr::slice_head(
    n = 20
  ) %>%
  
  print(
    n = 20,
    width = Inf
  )


# ============================================================
# STEP 3: REACTOME
# ============================================================

# Recreate Reactome collection if necessary

if (!exists("reactome_pathways")) {
  
  reactome_msig <- msigdbr::msigdbr(
    
    species = "Homo sapiens",
    
    collection = "C2",
    
    subcollection = "CP:REACTOME"
    
  )
  
  
  reactome_pathways <- split(
    
    reactome_msig$gene_symbol,
    
    reactome_msig$gs_name
    
  )
}


# ------------------------------------------------------------
# Run Reactome GSEA
# ------------------------------------------------------------

set.seed(12345)

fgsea_reactome <- fgsea::fgseaMultilevel(
  
  pathways = reactome_pathways,
  
  stats = gsea_ranks,
  
  minSize = 10,
  
  maxSize = 500,
  
  eps = 0
)


fgsea_reactome <- fgsea_reactome %>%
  
  tibble::as_tibble() %>%
  
  dplyr::mutate(
    
    direction =
      dplyr::if_else(
        NES > 0,
        "Older",
        "Younger"
      )
    
  ) %>%
  
  dplyr::arrange(
    padj
  )


# ------------------------------------------------------------
# Significant Reactome pathways
# ------------------------------------------------------------

reactome_sig <- fgsea_reactome %>%
  
  dplyr::filter(
    padj < 0.05
  ) %>%
  
  dplyr::arrange(
    padj
  )


# ------------------------------------------------------------
# Print Reactome summary
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "FINAL REACTOME GSEA SUMMARY\n",
  "============================================\n"
)

cat(
  "Pathways tested:",
  nrow(fgsea_reactome),
  "\n"
)

cat(
  "FDR < 0.05:",
  nrow(reactome_sig),
  "\n"
)

cat(
  "Older enriched:",
  sum(
    reactome_sig$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Younger enriched:",
  sum(
    reactome_sig$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nTop 30 Reactome pathways:\n"
)

fgsea_reactome %>%
  
  dplyr::select(
    pathway,
    size,
    NES,
    pval,
    padj,
    direction
  ) %>%
  
  dplyr::slice_head(
    n = 30
  ) %>%
  
  print(
    n = 30,
    width = Inf
  )


# ============================================================
# STEP 4: GO BIOLOGICAL PROCESS
# ============================================================

# Recreate GO BP collection if necessary

if (!exists("gobp_pathways")) {
  
  gobp_msig <- msigdbr::msigdbr(
    
    species = "Homo sapiens",
    
    collection = "C5",
    
    subcollection = "GO:BP"
    
  )
  
  
  gobp_pathways <- split(
    
    gobp_msig$gene_symbol,
    
    gobp_msig$gs_name
    
  )
}


# ------------------------------------------------------------
# Run GO BP GSEA
# ------------------------------------------------------------

set.seed(12345)

fgsea_gobp <- fgsea::fgseaMultilevel(
  
  pathways = gobp_pathways,
  
  stats = gsea_ranks,
  
  minSize = 10,
  
  maxSize = 500,
  
  eps = 0
)


fgsea_gobp <- fgsea_gobp %>%
  
  tibble::as_tibble() %>%
  
  dplyr::mutate(
    
    direction =
      dplyr::if_else(
        NES > 0,
        "Older",
        "Younger"
      )
    
  ) %>%
  
  dplyr::arrange(
    padj
  )


# ------------------------------------------------------------
# Significant GO BP pathways
# ------------------------------------------------------------

gobp_sig <- fgsea_gobp %>%
  
  dplyr::filter(
    padj < 0.05
  ) %>%
  
  dplyr::arrange(
    padj
  )


# ------------------------------------------------------------
# Print GO BP summary
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "FINAL GO BP GSEA SUMMARY\n",
  "============================================\n"
)

cat(
  "Pathways tested:",
  nrow(fgsea_gobp),
  "\n"
)

cat(
  "FDR < 0.05:",
  nrow(gobp_sig),
  "\n"
)

cat(
  "Older enriched:",
  sum(
    gobp_sig$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Younger enriched:",
  sum(
    gobp_sig$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nTop 40 GO BP pathways:\n"
)

fgsea_gobp %>%
  
  dplyr::select(
    pathway,
    size,
    NES,
    pval,
    padj,
    direction
  ) %>%
  
  dplyr::slice_head(
    n = 40
  ) %>%
  
  print(
    n = 40,
    width = Inf
  )


# ============================================================
# STEP 5: SEARCH OUR BIOLOGICAL THEMES
# ============================================================

pathway_keywords <- paste(
  
  c(
    "CYTOCHROME",
    "ARACHID",
    "EICOS",
    "PROSTAGLAND",
    "LEUKOTRIENE",
    "LIPOXYGEN",
    "FATTY_ACID",
    "LIPID",
    "ESTROGEN",
    "PPAR",
    "PEROXIS",
    "XENOBIOTIC",
    "OXIDATIVE",
    "REACTIVE_OXYGEN",
    "GLUTATHIONE",
    "MITOCHON",
    "CARDIAC",
    "MUSCLE"
  ),
  
  collapse = "|"
)


# ------------------------------------------------------------
# Hallmark theme hits
# ------------------------------------------------------------

hallmark_theme_hits <- fgsea_hallmark %>%
  
  dplyr::filter(
    
    grepl(
      pathway_keywords,
      pathway,
      ignore.case = TRUE
    )
    
  ) %>%
  
  dplyr::arrange(
    padj
  )


# ------------------------------------------------------------
# Reactome theme hits
# ------------------------------------------------------------

reactome_theme_hits <- fgsea_reactome %>%
  
  dplyr::filter(
    
    grepl(
      pathway_keywords,
      pathway,
      ignore.case = TRUE
    )
    
  ) %>%
  
  dplyr::arrange(
    padj
  )


# ------------------------------------------------------------
# GO BP theme hits
# ------------------------------------------------------------

gobp_theme_hits <- fgsea_gobp %>%
  
  dplyr::filter(
    
    grepl(
      pathway_keywords,
      pathway,
      ignore.case = TRUE
    )
    
  ) %>%
  
  dplyr::arrange(
    padj
  )


# ------------------------------------------------------------
# Print Hallmark theme hits
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "FINAL BIOLOGICALLY RELEVANT HALLMARK RESULTS\n",
  "============================================\n"
)

hallmark_theme_hits %>%
  
  dplyr::select(
    pathway,
    NES,
    pval,
    padj,
    direction
  ) %>%
  
  print(
    n = Inf,
    width = Inf
  )


# ------------------------------------------------------------
# Print Reactome theme hits
# ------------------------------------------------------------

cat(
  "\n============================================\n",
  "FINAL BIOLOGICALLY RELEVANT REACTOME RESULTS\n",
  "============================================\n"
)

reactome_theme_hits %>%
  
  dplyr::select(
    pathway,
    NES,
    pval,
    padj,
    direction
  ) %>%
  
  print(
    n = Inf,
    width = Inf
  )


# ============================================================
# STEP 6: CLEAN HALLMARK FIGURE
# ============================================================
#
# Use the top 12 pathways instead of 20 to make the
# presentation figure easier to read.
#
# ============================================================

hallmark_plot_data <- fgsea_hallmark %>%
  
  dplyr::filter(
    !is.na(padj)
  ) %>%
  
  dplyr::arrange(
    padj
  ) %>%
  
  dplyr::slice_head(
    n = 12
  ) %>%
  
  dplyr::mutate(
    
    pathway_clean =
      gsub(
        "^HALLMARK_",
        "",
        pathway
      ),
    
    pathway_clean =
      gsub(
        "_",
        " ",
        pathway_clean
      ),
    
    pathway_clean =
      factor(
        pathway_clean,
        levels = rev(pathway_clean)
      ),
    
    significance =
      dplyr::if_else(
        padj < 0.05,
        "FDR < 0.05",
        "FDR >= 0.05"
      )
  )


p_hallmark_gsea_final <- ggplot2::ggplot(
  
  hallmark_plot_data,
  
  ggplot2::aes(
    x = NES,
    y = pathway_clean
  )
  
) +
  
  ggplot2::geom_vline(
    
    xintercept = 0,
    
    linetype = "dashed",
    
    linewidth = 0.5
    
  ) +
  
  ggplot2::geom_point(
    
    ggplot2::aes(
      size = -log10(padj),
      shape = significance
    )
    
  ) +
  
  ggplot2::labs(
    
    title =
      "Hallmark pathways associated with age",
    
    subtitle =
      "Female LV ventricular cardiomyocytes; positive NES indicates enrichment in older females",
    
    x =
      "Normalized enrichment score (NES)",
    
    y =
      NULL,
    
    size =
      "-log10(FDR)",
    
    shape =
      "Evidence"
    
  ) +
  
  ggplot2::theme_classic(
    base_size = 13
  ) +
  
  ggplot2::theme(
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 16
      ),
    
    plot.subtitle =
      ggplot2::element_text(
        size = 12
      ),
    
    axis.text.y =
      ggplot2::element_text(
        size = 10
      ),
    
    legend.position =
      "bottom"
  )


print(
  p_hallmark_gsea_final
)


# ============================================================
# STEP 7: SAVE FINAL CORRECTED RESULTS
# ============================================================

dir.create(
  "results/enrichment",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Corrected ranked gene table
# ------------------------------------------------------------

readr::write_csv(
  
  gsea_rank_table,
  
  "results/enrichment/VCM_age_GSEA_ranked_genes_FINAL.csv"
  
)


# ------------------------------------------------------------
# Complete enrichment tables
# ------------------------------------------------------------

readr::write_csv(
  
  fgsea_hallmark,
  
  "results/enrichment/VCM_age_GSEA_Hallmark_all_FINAL.csv"
  
)


readr::write_csv(
  
  fgsea_reactome,
  
  "results/enrichment/VCM_age_GSEA_Reactome_all_FINAL.csv"
  
)


readr::write_csv(
  
  fgsea_gobp,
  
  "results/enrichment/VCM_age_GSEA_GO_BP_all_FINAL.csv"
  
)


# ------------------------------------------------------------
# Significant pathway tables
# ------------------------------------------------------------

readr::write_csv(
  
  hallmark_sig,
  
  "results/enrichment/VCM_age_GSEA_Hallmark_FDR05_FINAL.csv"
  
)


readr::write_csv(
  
  reactome_sig,
  
  "results/enrichment/VCM_age_GSEA_Reactome_FDR05_FINAL.csv"
  
)


readr::write_csv(
  
  gobp_sig,
  
  "results/enrichment/VCM_age_GSEA_GO_BP_FDR05_FINAL.csv"
  
)


# ------------------------------------------------------------
# Theme-specific tables
# ------------------------------------------------------------

readr::write_csv(
  
  hallmark_theme_hits,
  
  "results/enrichment/VCM_age_GSEA_Hallmark_target_themes_FINAL.csv"
  
)


readr::write_csv(
  
  reactome_theme_hits,
  
  "results/enrichment/VCM_age_GSEA_Reactome_target_themes_FINAL.csv"
  
)


readr::write_csv(
  
  gobp_theme_hits,
  
  "results/enrichment/VCM_age_GSEA_GO_BP_target_themes_FINAL.csv"
  
)


# ------------------------------------------------------------
# Save final clean Hallmark figure
# ------------------------------------------------------------

ggplot2::ggsave(
  
  filename =
    "results/figures/VCM_age_GSEA_Hallmark_FINAL.png",
  
  plot =
    p_hallmark_gsea_final,
  
  width =
    10,
  
  height =
    6.5,
  
  dpi =
    300
  
)


ggplot2::ggsave(
  
  filename =
    "results/figures/VCM_age_GSEA_Hallmark_FINAL.pdf",
  
  plot =
    p_hallmark_gsea_final,
  
  width =
    10,
  
  height =
    6.5
  
)


# ============================================================
# STEP 8: FINAL SUMMARY
# ============================================================

cat(
  "\n\n============================================\n",
  "FINAL CORRECTED GSEA COMPLETE\n",
  "============================================\n"
)


cat(
  "\nHALLMARK\n"
)

cat(
  "Significant pathways:",
  nrow(hallmark_sig),
  "\n"
)

cat(
  "Older:",
  sum(
    hallmark_sig$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Younger:",
  sum(
    hallmark_sig$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nREACTOME\n"
)

cat(
  "Significant pathways:",
  nrow(reactome_sig),
  "\n"
)

cat(
  "Older:",
  sum(
    reactome_sig$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Younger:",
  sum(
    reactome_sig$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nGO BIOLOGICAL PROCESS\n"
)

cat(
  "Significant pathways:",
  nrow(gobp_sig),
  "\n"
)

cat(
  "Older:",
  sum(
    gobp_sig$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Younger:",
  sum(
    gobp_sig$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\n============================================\n",
  "GSEA IS NOW FINAL\n",
  "============================================\n"
)