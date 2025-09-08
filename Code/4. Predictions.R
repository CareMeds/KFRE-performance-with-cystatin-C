################################################################################
## This script is intended to predict probabilities for each individual ########
## Author: Malou Magnani #######################################################
################################################################################
# remove history
rm(list=ls(all.names=TRUE))

# set seed for reproducibility
set.seed(27)

# set directory to load and save data
setwd("~/Data/")

# load data sets
load("cohort_outcomes.RData")

# load libraries
library(cmprsk)         # competing risk
library(dplyr)          # data manipulation

# load functions
source("~/Code/Functions for analyses.R")

################################################################################
### Describe outcome rates #####################################################
################################################################################
# Outcome categories: KF = 1, death = 2, censored = 0
cat("Number of individuals censored due to either emigration or end of follow up period within 2 years:", 
    sum(cohort$outcome_2y == 0), "\n",
    "Number of individuals censored due to emigration within 2 years:",
    sum(cohort$outcome_2y == 0 & 
          !is.na(cohort$date_emigration), na.rm = TRUE), "\n",
    "Number of individuals censored due to either emigration or end of follow up period within 5 years:", 
    sum(cohort$outcome_5y == 0), "\n",
    "Number of individuals censored due to emigration within 5 years:",
    sum(cohort$outcome_5y == 0 & 
          !is.na(cohort$date_emigration), na.rm = TRUE), "\n", 
    "\n",
    "Number of patients that have kidney failure within 2 years:",
    sum(cohort$outcome_2y == 1), "\n",
    "Number of patients that have kidney failure within 5 years:",
    sum(cohort$outcome_5y == 1), "\n", 
    "\n",
    "Number of patients that die within 2 years:",
    sum(cohort$outcome_2y == 2), "\n",
    "Number of patients that die within 5 years:",
    sum(cohort$outcome_5y == 2), "\n") 

################################################################################
### Cumulative incidence rates #################################################
################################################################################
cin_2y <- cmprsk::cuminc(ftime = cohort$time_to_event_2y, 
                         fstatus = cohort$outcome_2y, 
                         rho = 0, 
                         cencode = 0)
cin_5y <- cmprsk::cuminc(ftime = cohort$time_to_event_5y, 
                         fstatus = cohort$outcome_5y, 
                         rho = 0, 
                         cencode = 0)
cat(" Observed risk of KF at 2 years:             ", 
    cmprsk::timepoints(cin_2y, 2*365.25)$est[1,], "\n",
    "Observed risk of competing event at 2 years:", 
    cmprsk::timepoints(cin_2y, 2*365.25)$est[2,], "\n",
    "Observed risk of KF at 5 years:             ", 
    cmprsk::timepoints(cin_5y, 5*365.25)$est[1,], "\n",
    "Observed risk of competing event at 5 years:", 
    cmprsk::timepoints(cin_5y, 5*365.25)$est[2,], "\n")

################################################################################
### Modified time-to-event variable used for calculating modified C-statistic ##
################################################################################
# Wolbers et al.: to calculate the C-statistic whilst taking into account 
# competing risk, settime-to-event to infinity (i.e., prediction horizon) for 
# those experiencing the competing event, i.e. indicating that they will never 
# experience the outcome of interest (i.e., kidney failure). 
cohort$time_to_event_2y_modified <- ifelse(cohort$outcome_2y == 2, 
                                             2*365.25, 
                                             cohort$time_to_event_2y) 
cohort$time_to_event_5y_modified <- ifelse(cohort$outcome_5y == 2, 
                                             5*365.25, 
                                             cohort$time_to_event_5y) 

# censor competing event (i.e., death)
cohort <- cohort |> 
  dplyr::mutate(outcome_2y_modified = dplyr::case_when(outcome_2y == 0 ~ 0,
                                                       outcome_2y == 1 ~ 1,
                                                       outcome_2y == 2 ~ 0),
                outcome_5y_modified = dplyr::case_when(outcome_5y == 0 ~ 0,
                                                       outcome_5y == 1 ~ 1,
                                                       outcome_5y == 2 ~ 0)) 

################################################################################
### Prognostic index according to different equations ##########################
################################################################################
cohort <- cohort |>  
  dplyr::mutate(male = ifelse(female == 0, 1, 0), # create male variable
                alb = ifelse(alb == 0, 1, alb), # if alb = 0, change it to 1 (for log transformation)
                PI_ckd_epi_2009_cr = mapply(PI_KFRE, age, male, ckd_epi_2009_cr, alb), 
                PI_ckd_epi_2021_cr = mapply(PI_KFRE, age, male, ckd_epi_2021_cr, alb), 
                PI_ckd_epi_2012_cys = mapply(PI_KFRE, age, male, ckd_epi_2012_cys, alb),
                PI_ckd_epi_2012_cr_cys = mapply(PI_KFRE, age, male, ckd_epi_2012_cr_cys, alb),  
                PI_ckd_epi_2021_cr_cys = mapply(PI_KFRE, age, male, ckd_epi_2021_cr_cys, alb))

################################################################################
### Predicted 2-year KFRE risk according to different equations ################
################################################################################
cohort <- cohort |>
  dplyr::mutate(risk_2y_ckd_epi_2009_cr = mapply(KFRE_risk_2y, PI_ckd_epi_2009_cr), 
                risk_2y_ckd_epi_2021_cr = mapply(KFRE_risk_2y, PI_ckd_epi_2021_cr), 
                risk_2y_ckd_epi_2012_cys = mapply(KFRE_risk_2y, PI_ckd_epi_2012_cys),
                risk_2y_ckd_epi_2012_cr_cys = mapply(KFRE_risk_2y, PI_ckd_epi_2012_cr_cys), 
                risk_2y_ckd_epi_2021_cr_cys = mapply(KFRE_risk_2y, PI_ckd_epi_2021_cr_cys))

################################################################################
### Predicted 5-year KFRE risk according to different equations ################
################################################################################
cohort <- cohort |>  
  dplyr::mutate(risk_5y_ckd_epi_2009_cr = mapply(KFRE_risk_5y, PI_ckd_epi_2009_cr), 
                risk_5y_ckd_epi_2021_cr = mapply(KFRE_risk_5y, PI_ckd_epi_2021_cr), 
                risk_5y_ckd_epi_2012_cys = mapply(KFRE_risk_5y, PI_ckd_epi_2012_cys),
                risk_5y_ckd_epi_2012_cr_cys = mapply(KFRE_risk_5y, PI_ckd_epi_2012_cr_cys), 
                risk_5y_ckd_epi_2021_cr_cys = mapply(KFRE_risk_5y, PI_ckd_epi_2021_cr_cys))

# save cohort with predicted probabilities
save(cohort, cin_2y, cin_5y, file = "cohort_predictions.RData")