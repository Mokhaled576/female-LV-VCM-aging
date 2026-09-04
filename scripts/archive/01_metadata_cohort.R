# install.packages(c("tidyverse", "janitor", "here"))
install.packages(c("tidyverse", "janitor", "here"))
library(tidyverse)
library(janitor)
library(here)

donor_meta <- read_csv(
  "42003_2024_6582_MOESM18_ESM.csv",
  show_col_types = FALSE
)

dim(donor_meta)
names(donor_meta)
head(donor_meta)

# ------------------------------------------------------------
# Clean column names
# ------------------------------------------------------------

donor_meta <- donor_meta %>%
  clean_names()

names(donor_meta)

# ------------------------------------------------------------
# Define age groups
# ------------------------------------------------------------

donor_meta <- donor_meta %>%
  mutate(
    age_group = case_when(
      age < 20 ~ "Under_20",
      age >= 20 & age <= 45 ~ "Younger_20_45",
      age >= 46 & age <= 54 ~ "Transition_46_54",
      age >= 55 ~ "Older_55plus",
      TRUE ~ NA_character_
    )
  )

donor_meta %>%
  count(sex, age_group)

female_cohort <- donor_meta %>%
  filter(
    sex == "F",
    age_group %in% c("Younger_20_45", "Older_55plus")
  )

female_cohort %>%
  count(age_group)

female_cohort %>%
  arrange(age) %>%
  select(
    donor,
    age,
    age_group,
    data_source,
    sites_sampled
  )

sort(unique(donor_meta$sites_sampled))

female_cohort %>%
  count(data_source, age_group)

donor_meta %>%
  count(sex, age_group)

female_cohort %>%
  count(age_group)

sort(unique(donor_meta$sites_sampled))

female_cohort %>%
  count(data_source, age_group)