################################################################################
## This script is intended to create a flow chart of inclusion and exclusion ###
## Author: Malou Magnani #######################################################
################################################################################
# remove history
rm(list=ls(all.names=TRUE))

# set seed for reproducibility
set.seed(27)

# set directory to load and save data
setwd("~/Data/")

# load data
load("albuminuria_lab_test_clean.Rda")
load("death.Rda")
load("demographics.Rda")
load("lab_values.Rda")
load("rrt.Rda")

# load libraries
library(dplyr)          # data manipulation

# load functions
source("~/Code/Functions eGFR equations.R")

################################################################################
### Initial inclusion ##########################################################
################################################################################

# include patients with outpatient Cystatin C after 2011
cystatinc <- lab_values |> 
  dplyr::filter(test == "cystC") |>  
  dplyr::filter(ip == 0) |>
  dplyr::mutate(datum = as.Date(datum, format = "%Y-%m-%d")) |>  
  dplyr::select(lopnr, datum, result) |>
  dplyr::rename(cys = result, date_cys = datum) 
cystatinc_2011 <- cystatinc |> 
  dplyr::filter(date_cys >= "2011-01-01")

################################################################################
### Creatinine and Cystatin C measurements #####################################
################################################################################

# collect outpatient creatinine
creatinine <- lab_values |>  
  dplyr::filter(test == "crea") |> 
  dplyr::mutate(datum = as.Date(datum, format = "%Y-%m-%d")) |> 
  dplyr::rename(date_creat = datum) |>
  dplyr::rename(unit_creat = enhet) |>
  dplyr::rename(creat = result) |>
  dplyr::filter(ip==0) |> 
  dplyr::select(lopnr, date_creat, test,  unit_creat, creat) |> 
  dplyr::arrange(lopnr, date_creat)

# join the two data tables with creatinine and Cystatin C
cys_creat <- cystatinc_2011 |> 
  dplyr::left_join(creatinine) |> 
  dplyr::arrange(lopnr, date_creat, date_cys) |> 
  dplyr::filter(!is.na(creat)) # leave out missing creatinine

# keep only the rows if creatinine and Cystatin C were measured on the same day
cys_creat <- cys_creat |> 
  dplyr::filter(date_creat == date_cys) 

# count the number of measurements that were excluded 
# since creatinine and Cystatin C were not measured on the same day
removed_creat <- cystatinc_2011 |> 
  dplyr::filter(!(lopnr %in% cys_creat$lopnr))

################################################################################
### Select patients with Albuminuria measurement within 12 months of creatinine
### and Cystatin C measurement #################################################
################################################################################

# collect outpatient Albuminuria
alb <- albuminuria_lab_test_clean |>
  dplyr::filter(ip == 0) |>  
  dplyr::select(lopnr, date_alb = datum, alb = result, enhet)

# for each albuminuria measurement, 
# include only the closest creatinine/cystatin C measurement
alb_cys_creat <- merge(alb, cys_creat, by = "lopnr") |>  
  
  # calculate difference between albuminuria and creatinine date
  dplyr::mutate(diff = abs(as.numeric(as.Date(date_alb) - as.Date(date_creat)))) |>  
  
  # keep only row closer than 1 year
  dplyr::filter(diff < 365) |> 
  
  # for each measurement of albuminuria
  # keep the closest creatinine/Cystatin C measurement 
  dplyr::group_by(lopnr, date_alb) |> 
  dplyr::filter(diff == min(diff)) |> 
  
  #remove missing albuminuria
  dplyr::filter(!is.na(alb)) |> 
  
  dplyr::ungroup()

# count the number of patients that were removed
removed_alb <- cys_creat |> 
  dplyr::filter(!(lopnr %in% alb_cys_creat$lopnr))

################################################################################
### Remove registration anomalies ##############################################
################################################################################

# remove duplicate rows
death <- death |>
  dplyr::arrange(lopnr) |> 
  dplyr::distinct(lopnr, .keep_all = TRUE) |>
  dplyr::select(lopnr, death_date = dodsdat)

# remove measurements occurring after death
alb_cys_creat_death <- alb_cys_creat |> 
  dplyr::left_join(death, "lopnr") |> 
  dplyr::filter(is.na(death_date) | date_cys <= death_date & date_alb <= death_date)

# count the number of patients that were removed
removed_death <- alb_cys_creat |> 
  dplyr::filter(!(lopnr %in% alb_cys_creat_death$lopnr))

################################################################################
### Exclude age < 18 ###########################################################
################################################################################

# define baseline date as the first date at which albuminuria, creatinine, or 
# Cystatin C was measured
alb_cys_creat <- alb_cys_creat_death |> 
  dplyr::mutate(date_baseline = as.Date(ifelse(as.Date(date_alb) < as.Date(date_creat), 
                                        as.Date(date_creat), 
                                        as.Date(date_alb)))) 

# exclude those aged under 18
alb_cys_creat_18 <- alb_cys_creat |> 
  dplyr::left_join(demo,"lopnr") |> 
  dplyr::mutate(dob = as.Date(dob, format = "%Y-%m-%d")) |> 
  dplyr::mutate(age = round(as.numeric(date_creat-dob)/365.25,0)) |> 
  dplyr::select(-dob) |> 
  dplyr::filter(age>=18) |> 
  dplyr::arrange(lopnr)

# count number of patients that were under 18
removed_18 <- alb_cys_creat |> 
  dplyr::filter(!(lopnr %in% alb_cys_creat_18$lopnr))

################################################################################
### Exclude those eGFR < 10 and eGFR >= 60 #####################################
################################################################################

# exclude those eGFR < 10 and eGFR >= 60
alb_cys_creat_18_egfr <- alb_cys_creat_18 |> 
  dplyr::mutate(ckd_epi_2009_cr = ckd_epi_2009_cr(creat, age, female)) |> 
  dplyr::filter(ckd_epi_2009_cr < 60 & ckd_epi_2009_cr >= 10)

# count number of patients excluded
removed_egfr <- alb_cys_creat_18 |> 
  dplyr::filter(!(lopnr %in% alb_cys_creat_18_egfr$lopnr))

################################################################################
### Exclude those with RRT before baseline #####################################
################################################################################
cc_rrt <- alb_cys_creat_18_egfr |> 
  dplyr::left_join(rrt, "lopnr") |> 
  dplyr::mutate(new_rrt = ifelse(event_type == "GRAFT FAILURE"| 
                            event_type == "HD" | 
                            event_type == "PD" | 
                            event_type == "TX", 1, 0)) |>
  dplyr::filter((new_rrt == 0) |
           (new_rrt == 1 & date_baseline < rrt_date) |
           (is.na(rrt_date))) |> 
  dplyr::mutate(rrt_date = ifelse(new_rrt == 0, NA, rrt_date)) |> 
  dplyr::select(-event_type, -rrt)

# count number of patients with RRT before baseline
removed_krt <- alb_cys_creat_18_egfr |> 
  dplyr::filter(!(lopnr %in% cc_rrt$lopnr))

################################################################################
### Final flow chart ###########################################################
################################################################################
cat("Initial number of individuals included:", length(unique(cystatinc_2011$lopnr)), "\n",
    "Corresponding number of measurements:", length(cystatinc_2011$lopnr), "\n",
    "Patient exclusions due to registration anomalies:",
    length(unique(removed_death$lopnr)), "\n",
    "Corresponding number of measurements:", 
    length(removed_death$lopnr), "\n",
    "Patient exclusions due to criteria of having measured creatinine and Cystatin C on the same day:",
    length(unique(removed_creat$lopnr)), "\n",
    "Corresponding number of measurements:", length(removed_creat$lopnr), "\n",
    "Patient exclusions due to albuminuria measurement longer than 12 months:",
    length(unique(removed_alb$lopnr)), "\n",
    "Corresponding number of measurements:", 
    length(removed_alb$lopnr), "\n", 
    "Exclude patients under 18 years old:",
    length(unique(removed_18$lopnr)), "\n",
    "Corresponding number of measurements:",
    length(removed_18$lopnr), "\n",
    "Exclude patients with eGFR < 10 and EGFR >= 60:",
    length(unique(removed_egfr$lopnr)), "\n",
    "Corresponding number of measurements:",
    length(removed_egfr$lopnr), "\n",
    "Exclude patients with RRT before baseline:",
    length(unique(removed_krt$lopnr)), "\n",
    "Corresponding measurements:",
    length(removed_krt$lopnr), "\n",
    "Final number of inclusions:",
    length(unique(cc_rrt$lopnr)), "\n")

################################################################################
### Mean Albuminuria ###########################################################
################################################################################

# If multiple alb measurements on the same day
# take the mean measurement of albuminuria
data <- cc_rrt |>
  dplyr::group_by(lopnr, date_alb) |>  # operate per patient and date of the assay
  dplyr::mutate(alb = mean(alb)) |>    # if there are several values, take the mean on them
  dplyr::slice(1) |>                   # keep only one row per patient and date
  dplyr::ungroup()

################################################################################
### Choose one random measurement of each patient ##############################
################################################################################

# Choose one measurement date as index date per patient at random
cohort <- data |>
  # choose 1 random index date per patient
  dplyr::group_by(lopnr) |>
  dplyr::sample_n(1) |>
  dplyr::ungroup()

# exclude unnecessary variables
cohort <- cohort |> 
  dplyr::select(-dob_c, -dob_ym, -rrt_date_c, -diff, -test, -dob_c, -dob_ym)

################################################################################
### Event rate #################################################################
################################################################################
# 1626 KRT outcomes 
sum(cohort$new_rrt == 1, na.rm = TRUE) 

# save final cohort

save(cohort, file = "cohort.RData")
