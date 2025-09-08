################################################################################
## This script is intended to derive outcomes ##################################
## Author: Malou Magnani #######################################################
################################################################################
# remove history
rm(list=ls(all.names=TRUE))

# set seed for reproducibility
set.seed(27)

# set directory to load and save data
setwd("~/Data/")

# load data sets
load("cohort_covariates.RData")

# load libraries
library(dplyr)          # data manipulation

################################################################################
### Create time-to-event outcome for 2 and 5 years #############################
################################################################################
cohort <- cohort |> 
  dplyr::mutate(
    # define two years in the future from index date
    end_follow_up_2y = as.Date(index_dt + 365.25 * 2),
    
    # define date at which event occurs 
    date_2y = pmin(as.Date(rrt_date),         # outcome of interest
                   as.Date(death_date),       # competing event
                   as.Date("2021-12-31"),     # censoring
                   as.Date(date_emigration),  # censoring 
                   as.Date(end_follow_up_2y), na.rm = TRUE), # censoring
    
    # create a time-to-event variable
    time_to_event_2y =  as.numeric(date_2y - index_dt),
    
    # create indicator variable for event
    outcome_2y = dplyr::case_when(date_2y == rrt_date ~ 1,   # outcome of interest
                           date_2y == death_date ~ 2, # competing event
                           TRUE ~ 0),              # censoring
    
    # define five years in the future from index date
    end_follow_up_5y = as.Date(index_dt + 365.25 * 5),
    
    # define date at which event occurs
    date_5y = pmin(as.Date(rrt_date), 
                   as.Date(death_date), 
                   as.Date("2021-12-31"), 
                   as.Date(date_emigration), 
                   as.Date(end_follow_up_5y), na.rm = TRUE),
    
    # create a time-to-event variable
    time_to_event_5y =  as.numeric(date_5y - index_dt),
    
    # create indicator variable for event
    outcome_5y = dplyr::case_when(date_5y == rrt_date ~ 1,   # outcome of interest
                           date_5y == death_date ~ 2, # competing event
                           TRUE ~ 0))                 # censoring

# save cohort with outcomes
save(cohort, file = "cohort_outcomes.RData")