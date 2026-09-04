# ============================================================
# 01_metadata_cohort.R
# Female non-failing LV donor cohort
# ============================================================
#
# Source:
# Read et al. (2024), Communications Biology
# DOI: 10.1038/s42003-024-06582-y
#
# Primary comparison:
#   Younger: 20-45 years
#   Older:   >=55 years
#
# Exclusions:
#   <20 years
#   46-54 years
#
# IMPORTANT:
# Older_55plus is an age-based proxy and is NOT confirmed
# menopausal status.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

donor_file <- file.path(
  "data_raw",
  "42003_2024_6582_MOESM18_ESM.csv"
)

stopifnot(file.exists(donor_file))

donor_meta <- readr::read_csv(
  donor_file,
  show_col_types = FALSE
) %>%
  janitor::clean_names()

required_cols <- c(
  "donor", "age", "sex", "data_source", "sites_sampled"
)

stopifnot(all(required_cols %in% names(donor_meta)))

# ------------------------------------------------------------
# Define age groups
# ------------------------------------------------------------

donor_meta <- donor_meta %>%
  dplyr::mutate(
    age_group = dplyr::case_when(
      age < 20 ~ "Under_20",
      age >= 20 & age <= 45 ~ "Younger_20_45",
      age >= 46 & age <= 54 ~ "Transition_46_54",
      age >= 55 ~ "Older_55plus",
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------
# Female donors in the two prespecified age groups
# ------------------------------------------------------------

female_age_cohort <- donor_meta %>%
  dplyr::filter(
    sex == "F",
    age_group %in% c("Younger_20_45", "Older_55plus")
  )

# ------------------------------------------------------------
# Primary tissue: left ventricle
#
# Sites_Sampled may contain comma-separated sites such as
# "Apex,LV". Match LV as a complete comma-delimited token.
# ------------------------------------------------------------

female_lv <- female_age_cohort %>%
  dplyr::filter(
    stringr::str_detect(
      sites_sampled,
      stringr::regex("(^|,)\\s*LV\\s*(,|$)", ignore_case = TRUE)
    )
  ) %>%
  dplyr::arrange(age)

# ------------------------------------------------------------
# Cohort checks
# ------------------------------------------------------------

cat("\nFemale LV donor cohort\n")
cat("----------------------\n")
cat("Total donors:", nrow(female_lv), "\n")

print(
  female_lv %>%
    dplyr::count(age_group, name = "n_donors")
)

print(
  female_lv %>%
    dplyr::count(data_source, age_group, name = "n_donors")
)

stopifnot(sum(duplicated(female_lv$donor)) == 0)
stopifnot(sum(is.na(female_lv$donor)) == 0)
stopifnot(sum(is.na(female_lv$age_group)) == 0)

# Expected for the primary analysis:
#   24 donors total
#   5 younger
#   19 older

# ------------------------------------------------------------
# Save cohort tables
# ------------------------------------------------------------

readr::write_csv(
  female_lv,
  "results/tables/female_LV_primary_cohort.csv"
)

readr::write_lines(
  female_lv$donor,
  "results/tables/female_LV_donor_IDs.txt"
)

study_age_table <- female_lv %>%
  dplyr::count(data_source, age_group, name = "n_donors") %>%
  tidyr::pivot_wider(
    names_from = age_group,
    values_from = n_donors,
    values_fill = 0
  )

readr::write_csv(
  study_age_table,
  "results/tables/female_LV_study_age_counts.csv"
)

cat("\nSaved cohort outputs to results/tables/\n")
