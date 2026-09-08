################################################################################
## This script is intended to derive outcomes ##################################
## Author: Malou Magnani #######################################################
################################################################################
# remove history
rm(list=ls(all.names=TRUE))

# set seed for reproducibility
set.seed(27)

# set directory to load and save data
setwd("P:/SCREAM2/SCREAM2_Research/Malou Magnani/Final/Data/")

# load data sets
load("cohort_covariates.RData")

# load libraries
library(dplyr)          # data manipulation

################################################################################
### Create time-to-event outcomes for each horizon ############################
################################################################################
# named list of horizons (in years) to derive outcomes for. Inf = no horizon
# cap at all -- follow-up runs to the earliest of event/death/admin-censoring/
# emigration only. Useful for time-dependent AUC at a horizon t0 close to a
# capped end_follow_up, where truncating time-to-event exactly at t0 can tie
# the risk set and make it impossible to have any subject with T > t0.
horizons <- list("2y" = 2, "5y" = 5, "inf" = Inf)

for (suffix in names(horizons)) {
  horizon_years <- horizons[[suffix]]
  
  # date at which event occurs: earliest of outcome of interest, competing
  # event, or the censoring dates -- with the horizon cap included only when
  # horizon_years is finite
  censoring_dates <- list(
    as.Date(cohort$rrt_date),         # outcome of interest
    as.Date(cohort$death_date),       # competing event
    as.Date("2021-12-31"),            # censoring
    as.Date(cohort$date_emigration)   # censoring
  )
  if (is.finite(horizon_years)) {
    end_follow_up <- as.Date(cohort$index_dt + 365.25 * horizon_years)
    cohort[[paste0("end_follow_up_", suffix)]] <- end_follow_up
    censoring_dates <- c(censoring_dates, list(end_follow_up)) # censoring
  }
  date_event <- do.call(pmin, c(censoring_dates, na.rm = TRUE))
  
  # create a time-to-event variable
  cohort[[paste0("time_to_event_", suffix)]] <- as.numeric(date_event - cohort$index_dt)
  
  # create indicator variable for event
  cohort[[paste0("outcome_", suffix)]] <- dplyr::case_when(
    date_event == cohort$rrt_date ~ 1,   # outcome of interest
    date_event == cohort$death_date ~ 2, # competing event
    TRUE ~ 0                             # censoring
  )
}

# save cohort with outcomes
save(cohort, file = "cohort_outcomes.RData")