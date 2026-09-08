################################################################################
## This script is intended to validate predictions #############################
## Author: Malou Magnani and Carolien C.H.M. Maas ##############################
##
## This script runts:
### main results table -> calibration plots -> discrimination forest plot ->
## DCA -> risk distributions -> risk differences -> eGFR distributions ->
## reclassification tables -> sensitivity analyses (outcome comparison,
## eGFR subgroups, eGFR range)
################################################################################

# remove history
rm(list = ls(all.names = TRUE))

# set seed for reproducibility
set.seed(27)

# set directory to save results
setwd("P:/SCREAM2/SCREAM2_Research/Malou Magnani/Final/")

# load data sets
load("Data/cohort_predictions.RData")

# load functions
source("Code/Functions for analyses.R")

# load libraries
library(survival)  # time-to-event analyses
library(patchwork) # combine plots

################################################################################
### Output directory ###########################################################
################################################################################
# every table/plot from this script is written here, not into the working
# directory directly
data_dir <- "Data"
results_dir <- "Results"
dir.create(results_dir, showWarnings = FALSE)

################################################################################
### Shared constants ###########################################################
################################################################################
horizons <- c(2, 5)

equations <- c(
  "ckd_epi_2009_cr",
  "ckd_epi_2021_cr",
  "ckd_epi_2012_cys",
  "ckd_epi_2012_cr_cys",
  "ckd_epi_2021_cr_cys"
)

model_names <- c(
  "CKD-EPIcr 2009",
  "CKD-EPIcr 2021",
  "CKD-EPIcys 2012",
  "CKD-EPIcrcys 2012",
  "CKD-EPIcrcys 2021"
)

colors_palette <- c("darkorange4",
                    "darkred",
                    "darkorchid4",
                    "darkblue",
                    "darkgreen")

n_bootstraps <- 500

# 95th percentile of 5-year predicted risk (CKD-EPIcr 2009), used to set
# shared, "trimmed" axis limits across several plots.
max_5y <- quantile(cohort$risk_5y_ckd_epi_2009_cr, 0.95, na.rm = TRUE)

################################################################################
### 1. Main performance measures table -> Measures.xlsx ########################
################################################################################
## `compute_performance_measures()` now generates its own column layout
## internally (via `measure_colnames()`, shared with the sensitivity
## tables), so there's no `measures` argument to keep in sync by hand
## anymore -- see Functions_for_analyses.R for details.
measures_df <- compute_performance_measures(
  cohort = cohort,
  horizons = horizons,
  equations = equations,
  model_names = model_names,
  B = n_bootstraps,
  file = file.path(results_dir, "All_measures.xlsx")
)

################################################################################
### 2. Calibration plots (+ predicted-risk histogram insets) ###################
################################################################################
plot_calibration_curves(
  cohort = cohort,
  horizons = horizons,
  equations = equations,
  model_names = model_names,
  colors_palette = colors_palette,
  max_5y = max_5y,
  file_full = file.path(results_dir, "Combined calibration plot.png"),
  file_trimmed = file.path(results_dir, "Combined calibration plot trimmed.png")
)

################################################################################
### 3. Discrimination forest plot #################################£############
################################################################################
# measures_file must match the `file` used in step 1 above, since this reads
# the AUC + CI columns straight out of that Excel file
plot_discrimination_forest(
  measures_file = file.path(results_dir, "All_measures.xlsx"),
  out_file = file.path(results_dir, "forest_plot_AUC.png")
)

################################################################################
### 4. Decision curve analysis (DCA) ###########################################
################################################################################
# The raw DCA computation is slow, so it's cached; set new_DCA_run = TRUE to
# regenerate it (e.g. after the underlying predictions change).
run_dca(
  cohort = cohort,
  horizons = horizons,
  equations = equations,
  model_names = model_names,
  colors_palette = colors_palette,
  max_5y = max_5y,
  new_DCA_run = FALSE,
  cache_file = file.path(data_dir, "analyses/dca.Rda"),
  plot_file = file.path(results_dir, "Combined DCA plot.png"),
  table_file = file.path(results_dir, "Net Benefit Table.xlsx")
)

################################################################################
### 5. Predicted risk distributions ############################################
################################################################################
risk_2y_long <- build_risk_distribution_long(
  cohort,
  horizon = 2,
  equations = equations,
  model_names = model_names
)
risk_5y_long <- build_risk_distribution_long(
  cohort,
  horizon = 5,
  equations = equations,
  model_names = model_names
)

plot_risk_distributions(
  risk_2y_long,
  risk_5y_long,
  out_file = file.path(results_dir, "Combined distribution plots.png")
)

################################################################################
### 6. Risk differences vs. reference equation #################################
################################################################################
plot_risk_differences(
  risk_2y_long,
  risk_5y_long,
  reference_label = "CKD-EPIcr 2009",
  out_file = file.path(results_dir, "Combined risk difference.png")
)

################################################################################
### 7. eGFR distributions ######################################################
################################################################################
plot_egfr_distributions(
  cohort,
  equations = equations,
  model_names = model_names,
  out_file = file.path(results_dir, "eGFR distribution.png")
)
write_egfr_summary_table(
  cohort = cohort,
  equations = equations,
  out_file = file.path(results_dir, "eGFR distributions.xlsx")
)

################################################################################
### 8. Reclassification tables #################################################
################################################################################
build_reclassification_tables(
  cohort,
  equations = equations,
  file_2y = file.path(results_dir, "Reclassification 2-year predictions.xlsx"),
  file_5y = file.path(results_dir, "Reclassification 5-year predictions.xlsx")
)

################################################################################
### 9. Sensitivity analysis: outcome comparison ################################
################################################################################
run_sensitivity_outcome_comparison(
  cohort = cohort,
  horizons = horizons,
  equations = equations,
  model_names = model_names,
  out_file = file.path(results_dir, "SA_tdAUC_outcomes.xlsx"),
  B = n_bootstraps
)

################################################################################
### 10. Sensitivity analysis: eGFR subgroups ###################################
################################################################################
run_sensitivity_egfr_subgroup(
  cohort = cohort,
  horizons = horizons,
  equations = equations,
  model_names = model_names,
  out_file = file.path(results_dir, "SA_eGFR_subgroup.xlsx"),
  B = n_bootstraps
)

################################################################################
### 11. Sensitivity analysis: eGFR range (10-60) ###############################
################################################################################
run_sensitivity_egfr_range(
  cohort = cohort,
  horizons = horizons,
  equations = equations,
  model_names = model_names,
  out_file = file.path(results_dir, "SA_eGFR_range.xlsx"),
  B = n_bootstraps
)