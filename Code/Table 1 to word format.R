# Load required libraries
library(gtsummary)
library(dplyr)
library(flextable)
library(officer)
library(papaja)
library(tableone)
library(knitr)
library(tidyverse)
library(kableExtra)
library(gt)

# Load datasets
load("P:/SCREAM2/SCREAM2_Research/Malou Magnani/Final/Data/cohort_outcomes.Rdata")

# Combine datasets
cohort_outcomes_c <- cohort |>
  ungroup() |>
  mutate(creat_US = creat / 88.4)

# Combine the datasets, ensuring no column duplication
cohort_outcomes_c <- cohort_outcomes_c |> 
  mutate(
    egfr_category = case_when(
      ckd_epi_2009_cr >= 45 ~ "45 to <60",
      ckd_epi_2009_cr >= 30 & ckd_epi_2009_cr < 45 ~ "30 to <45",
      ckd_epi_2009_cr >= 15 & ckd_epi_2009_cr < 30 ~ "15 to <30",
      ckd_epi_2009_cr < 15 ~ "<15",
      TRUE ~ NA_character_
    ),
    egfr_category = factor(egfr_category, levels = c("45 to <60", "30 to <45", "15 to <30", "<15"))
  )


# Create the empty variables for table layout
cohort_outcomes_c$empty1 <- NA
cohort_outcomes_c$empty2 <- NA
cohort_outcomes_c$empty3 <- NA


# List of variables to include in the table
listvar <- c("age", "female", "egfr_category", "empty1", "ckd_epi_2009_cr", "ckd_epi_2021_cr", "ckd_epi_2012_cys", "ckd_epi_2012_cr_cys", "ckd_epi_2021_cr_cys", "creat", "cys", "alb", "empty2", "mi", "ihd", "hyperten", "hf", "stroke", "cevd", "arrh", "pvd", "dm", "cancer", "copd", "liver", "empty3", "bblock", "hypoglycemic", "ccb", "diur", "rasi", "lipid", "nsaid")

# Continuous and categorical variables
continuous <- c("age", "ckd_epi_2009_cr", "ckd_epi_2021_cr", "ckd_epi_2012_cys", "ckd_epi_2012_cr_cys", "ckd_epi_2021_cr_cys", "creat", "cys", "alb")
catvar <- listvar[!listvar %in% continuous]

# Define labels for the variables
labels <- c("Age, median (IQR), y", 
            "Female, n (%)", 
            "eGFR category, n (%)",
            "45 to <60",
            "30 to <45",
            "15 to <30",
            "10 to <15",
            "Median eGFR (IQR), mL/min/1.73m2",
            "  CKD-EPIcr 2009, ml/min/1.73m2",
            "  CKD-EPIcr 2021, ml/min/1.73m2",
            "  CKD-EPIcys 2012, ml/min/1.73m2",
            "  CKD-EPIcr-cys 2012, ml/min/1.73m2",
            "  CKD-EPIcr-cys 2021, ml/min/1.73m2",
            "Median serum creatinine (IQR), mg/dL",
            "Median serum cystatin C (IQR), mg/L",
            "Median UACR (IQR), mg/g",
            "Comorbidities, n (%)",
            " Myocardial infarction",
            " Other Ischemic Heart Disease",
            " Hypertension",
            " Heart Failure",
            " Stroke",
            " Other cerebrovascular disease",
            " Arrhythmia",
            " Peripheral vascular disease",
            " Diabetes mellitus",
            " Cancer in previous year",
            " Chronic obstructive pulmonary disease",
            " Liver disease",
            "Medications, n (%)",
            " Beta blocker",
            " Calcium channel blocker",
            " Diabetes medications",
            " Diuretic",
            " ACEi/ARB",
            " Lipid lowering drug",
            " NSAID")

# Construct mutually non-exclusive outcome groups
kfrt_2y <- cohort_outcomes_c |> 
  filter(outcome_2y == 1) |> 
  mutate(outcome_group = "KFRT within 2 years")

kfrt_5y <- cohort_outcomes_c |> 
  filter(outcome_5y == 1) |> 
  mutate(outcome_group = "KFRT within 5 years")

no_kfrt <- cohort_outcomes_c |> 
  filter(outcome_5y == 0 | outcome_5y == 2) |> 
  mutate(outcome_group = "No KFRT")

# Bind rows: allows patients to appear in multiple groups
combined_data <- bind_rows(kfrt_2y, kfrt_5y, no_kfrt)

# Create a deduplicated dataset for the "overall" group (1 row per patient)
overall_data <- cohort_outcomes_c |>
  distinct(lopnr, .keep_all = TRUE) |> 
  mutate(outcome_group = "Overall")

# Combine overall + grouped data
combined_with_overall <- bind_rows(overall_data, combined_data)

# Create the table
table_ovr_gtsummary <- combined_with_overall |>
  select(outcome_group, all_of(listvar)) |>
  tbl_summary(
    by = outcome_group,
    statistic = list(
      all_of(continuous) ~ "{median} ({p25}, {p75})",
      all_of(catvar) ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 0,
      all_categorical() ~ c(0, 0),
      "cys" ~ 2
    ),
    missing = "no"
  ) |>
  modify_table_body(~ .x |> mutate(across(everything(), ~ ifelse(. %in% "0 (NA%)", "", .)))) |>
  modify_table_body(~ .x |> mutate(label = labels)) |>
  modify_footnote(everything() ~ NA)


# Convert to gt table
table_ovr_gt <- as_gt(table_ovr_gtsummary) |>
  tab_header(title = "Baseline characteristics split by outcome")

# Save the table to a Word document
gtsave(table_ovr_gt, "P:/SCREAM2/SCREAM2_Research/Malou Magnani/Final/Results/baseline_characteristics_table.docx")