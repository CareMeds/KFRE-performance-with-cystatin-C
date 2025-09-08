################################################################################
## This script is intended to create a flow chart of inclusion and exclusion ###
## Author: Malou Magnani #######################################################
################################################################################
# remove history
rm(list=ls(all.names=TRUE))

# Load libraries
library(dplyr)      # data manipulation
library(gtsummary)  # summarize data
library(tidyselect) # select variables
library(gt)         # table header and save

# Load datasets
load("~/Data/cohort_predictions.Rdata")

# categorical eGFR
cohort <- cohort |>
  dplyr::mutate(cat_eGFR_high = ifelse(ckd_epi_2009_cr >= 45 & ckd_epi_2009_cr < 60, 1, 0),
                cat_eGFR_medh = ifelse(ckd_epi_2009_cr >= 30 & ckd_epi_2009_cr < 45, 1, 0),
                cat_eGFR_medl = ifelse(ckd_epi_2009_cr >= 15 & ckd_epi_2009_cr < 30, 1, 0),
                cat_eGFR_low = ifelse(ckd_epi_2009_cr >= 10 & ckd_epi_2009_cr < 15, 1, 0))

# Create the empty variables for table layout
cohort$empty1 <- NA
cohort$empty2 <- NA
cohort$empty3 <- NA
cohort$empty4 <- NA
cohort$empty5 <- NA
cohort$empty6 <- NA
cohort$empty7 <- NA

# List of variables to include in the table
listvar <- c("age", 
             "female", 
             "empty1", # eGFR categorical
             "cat_eGFR_high",
             "cat_eGFR_medh",
             "cat_eGFR_medl",
             "cat_eGFR_low",
             "empty2", # equations
             "ckd_epi_2009_cr", 
             "ckd_epi_2021_cr", 
             "ckd_epi_2012_cys", 
             "ckd_epi_2012_cr_cys", 
             "ckd_epi_2021_cr_cys", 
             "creat", "cys", "alb",
             "empty3", # comorbidities
             "mi", "ihd", "hyperten", "hf", "stroke", "cevd", "arrh", "pvd", 
             "dm", "cancer", "copd", "liver", 
             "empty4", # medications
             "bblock", "hypoglycemic", "ccb", "diur", "rasi", "lipid", "nsaid",
             "outcome_2y", "time_to_event_2y",
             "outcome_5y", "time_to_event_5y")

# Continuous and categorical variables
continuous <- c("age", "ckd_epi_2009_cr", "ckd_epi_2021_cr", 
                "ckd_epi_2012_cys", "ckd_epi_2012_cr_cys", "ckd_epi_2021_cr_cys", 
                "creat", "cys", "alb", "time_to_event_2y", "time_to_event_5y")
catvar <- listvar[!listvar %in% continuous]

# Define labels for the variables
labels <- c("Age, median (IQR), y", 
            "Female, n (%)", 
            "eGFR category, n(%)",
            "  45 to <60",
            "  30 to <45",
            "  15 to <30",
            "  10 to <15",
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
            " NSAID",
            "Outcome 2 years, n (%)",
            " Alive without KFRT",
            " KFRT",
            " Death without KFRT",
            " Median time to event (IQR), days",
            "Outcome 5 years, n (%)",
            " Alive without KFRT",
            " KFRT",
            " Death without KFRT",
            " Median time to event (IQR), days")

# Create the summary table using gtsummary
table_ovr_gtsummary <- cohort |>
  dplyr::ungroup() |>
  dplyr::select(tidyselect::all_of(listvar)) |>
  gtsummary::tbl_summary(
    by = NULL,  # No grouping
    statistic = list(
      tidyselect::all_of(continuous) ~ "{median} ({p25}, {p75})",  # Median and IQR for continuous
      tidyselect::all_of(catvar) ~ "{n} ({p}%)"                    # Count and percentage for categorical
    ),
    digits = list(
      gtsummary::all_continuous() ~ 0,         # Round continuous to 1 decimal
      gtsummary::all_categorical() ~ c(0, 1),  # Round percentages to 1 decimal
      "cys" ~ 2,
      "time_to_event_2y" ~ 1,
      "time_to_event_5y" ~ 1
    ),
    missing = "no") |>  # Exclude missing values
  gtsummary::modify_table_body(~ .x |> 
                                 dplyr::mutate(across(everything(), ~ ifelse(. %in% "0 (NA%)", "", .)))  # Replace values with empty strings
  ) |>
  gtsummary::modify_table_body(~ .x |> 
                                 dplyr::mutate(label = labels)) |> # Assign labels
  gtsummary::modify_footnote(everything() ~ NA)

# Convert to gt table
table_ovr_gt <- gtsummary::as_gt(table_ovr_gtsummary) |>
  gt::tab_header(title = "Table 1. Baseline characteristics")

# Save the table to a Word document
gt::gtsave(table_ovr_gt, "~/Results/Table 1. Baseline Characteristics.docx")