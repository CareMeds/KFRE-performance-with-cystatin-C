################################################################################
## This script is intended to obtain covariates ################################
## Author: Malou Magnani #######################################################
################################################################################
# remove history
rm(list=ls(all.names=TRUE))

# set seed for reproducibility
set.seed(27)

# set directory to load and save data
setwd("~/Data/")

# Load data sets
load("cohort.Rdata")
load("emigration.Rda")    # emigration
load("diagnoses_kon.Rda") # primary care
load("diagnoses_ovr.Rda") # specialist care
load("diagnoses_slv.Rda") # hospitalizations
load("medication.Rda")    # medications

# load libraries
library(dplyr)          # data manipulation

# load functions
source("~/Code/Functions eGFR equations.R")

################################################################################
### Obtain emigration data #####################################################
################################################################################
emigration <- emigration |> 
  dplyr::filter(hkod == "U") |> # select only emigration
  dplyr::mutate(hdat = as.Date(hdat, format = "%Y-%m-%d")) |> # format date
  dplyr::select(lopnr, hdat) |> # keep only patient ID (lopnr) and emigration date (hdat)
  dplyr::arrange(lopnr, hdat)

# Obtain the first date at which individuals emigrate
# Only keep the first emigration date after baseline
cohort <- cohort |> 
  dplyr::left_join(emigration, "lopnr") |> # add emigration to cohort
  dplyr::rename(date_emigration = hdat) |> # rename the date of emigration column
  # calculate difference between creat/index date and emigration date
  # if the emigration date is before the index date, emigration date is set to NA
  dplyr::mutate(diff = as.numeric(as.Date(date_emigration) - as.Date(date_baseline)), 
                date_emigration = ifelse(diff > 0, date_emigration, NA)) |> 
  dplyr:: arrange(lopnr, date_emigration) |> # earliest emigration dates (after index date) first
  dplyr::group_by(lopnr) |> 
  dplyr::slice(1) |> # select only the first emigration date after the index date per patient
  dplyr::mutate(date_emigration = as.Date(date_emigration))

################################################################################
### Calculate eGFR using different equations ###################################
################################################################################
cohort_egfr <- cohort |> 
  dplyr::mutate(ckd_epi_2009_cr = ckd_epi_2009_cr(creat, age, female),
                ckd_epi_2021_cr = ckd_epi_2021_cr(creat, age, female),
                ckd_epi_2012_cys = ckd_epi_2012_cys(cys, age, female),
                ckd_epi_2012_cr_cys = ckd_epi_2012_cr_cys(creat, cys, age, female),
                ckd_epi_2021_cr_cys = ckd_epi_2021_cr_cys(creat, cys, age, female))

# only keep necessary information
cohort <- cohort_egfr |>
  dplyr::arrange(lopnr, date_creat, date_cys, date_alb) |>
  dplyr::rename(index_dt = date_baseline) |>
  dplyr::select(lopnr, index_dt, female, age, creat, cys, alb, new_rrt, rrt_date, 
                death_date, date_emigration, ckd_epi_2009_cr, ckd_epi_2021_cr, 
                ckd_epi_2012_cys, ckd_epi_2012_cr_cys, ckd_epi_2021_cr_cys) 

################################################################################
### Exctact comorbidities diagnoses ############################################
################################################################################
# ICD code and date for primary care visit
diagnoses_kon <- diagnoses_kon |> 
  dplyr::filter(LOPNR %in% cohort$lopnr) |> # select patients included in current cohort
  dplyr::mutate(bdat = as.Date(bdat, format ="%Y-%m-%d")) |> 
  dplyr::rename(date_diag = bdat, lopnr = LOPNR) |> # date when ICD code was assigned 
  dplyr::select(lopnr, date_diag, diagnosis)

# ICD code and date for outpatient care
diagnoses_ovr <- diagnoses_ovr |> 
  dplyr::filter(lopnr %in% cohort$lopnr) |> 
  dplyr::mutate(bdat = as.Date(bdat, format ="%Y-%m-%d")) |> 
  dplyr::rename(date_diag = bdat) |> 
  dplyr::select(lopnr, date_diag, diagnosis)

# ICD code and date for hospital care
diagnoses_slv <- diagnoses_slv |> 
  dplyr::filter(lopnr %in% cohort$lopnr) |> 
  dplyr::mutate(bdat = as.Date(bdat, format ="%Y-%m-%d")) |> 
  dplyr::rename(date_diag = bdat) |> 
  dplyr::select(lopnr, date_diag, diagnosis)

# combine data sets from primary, outpatient, and hospital care
diagnoses <- rbind(diagnoses_kon,diagnoses_ovr,diagnoses_slv) |> 
  dplyr::arrange(lopnr, date_diag) |>  
  dplyr::select(lopnr, diagnosis, date_diag)

# obtain index_dt from cohort
diagnoses_index_dt <- diagnoses |>
  dplyr::left_join(cohort, by="lopnr") |>
  dplyr::select(lopnr, diagnosis, date_diag, index_dt)

save(diagnoses_index_dt, file = "diagnoses.Rdata")

################################################################################
### Add comorbidities any time before the index date to cohort #################
################################################################################
# select diagnosis per ID, ordered on date
diagnoses_per_ID <- diagnoses_index_dt |>
  dplyr::filter(date_diag <= index_dt &  # diagnosis before any time before baseline
                  !is.na(date_diag)) |>
  dplyr::group_by(lopnr, diagnosis) |>
  dplyr::arrange(date_diag) |>
  dplyr::slice(1) |>
  dplyr::select(lopnr, diagnosis) |>
  dplyr::rename(comorb = diagnosis)

# enlist the comorbidities ICD codes 
comorb_index <- data.frame(names=c("mi", "ihd", "hyperten", "hf", "stroke", "cevd", 
                                   "arrh", "pvd", "dm", "copd", "liver"),
                           ICD=c("^I200|^I21|^I22",            # myocaridal infarction
                                 "^I201|^I208|^I209|^I24|^I25",  # ischemic heart disease
                                 "^I10|^I11|^I12|^I13|^I14|^I15", # hypertension
                                 "^I110|^I130|^I132|^I50",      # heart failure
                                 "^I60|^I61|^I62|^I63|^I64|^I693|^I698|^I694", # stroke
                                 "^I65|^I66|^I67|^I68|^I69|^G450|^G451|^G452|^G453|^G458|^G459|^G46", # cerebrovascular disease
                                 "^I44|^I45|^I46|^I47|^I48|^I49", # arrhythmia 
                                 "^I70|^I72|^I73",             # peripheral vascular disease
                                 "^E10|^E11|^E12|^E13|^E14",     # diabetes mellitus
                                 "^J44",                     # chronic obstructive pulmonary disease
                                 "^B18|^I982|^K70|^K71|^K72|^K73|^K74|^K75|^K76|^K77")) # liver disease
for (ICD_index in 1:length(comorb_index$names)){
  # indicate for each ID if the ICD was present
  diagnosis_per_ICD <- diagnoses_per_ID |>
    dplyr::mutate(ICD_present = grepl(comorb_index$ICD[ICD_index], comorb)) |>
    dplyr::group_by(lopnr) |>
    dplyr::summarise(ICD_present = as.integer(any(ICD_present))) |> # keep one row per lopnr
    dplyr::rename(!!comorb_index$names[ICD_index] := ICD_present)   # rename column to current ICD
  
  # join with cohort
  cohort <- cohort |>
    dplyr::left_join(diagnosis_per_ICD, by="lopnr") |>
    dplyr::mutate(!!comorb_index$names[ICD_index] := # those who are not in diagnoses_per_ICD did not have that comorbidity
                    tidyr::replace_na(!!dplyr::sym(comorb_index$names[ICD_index]), 0))
}

################################################################################
### Add cancer diagnosis one year before the index date to cohort ##############
################################################################################
# select diagnosis per ID, ordered on date
diagnoses_per_ID_1y <- diagnoses_index_dt |>
  dplyr::filter(date_diag <= index_dt &         # diagnosis before baseline
                  date_diag >= (index_dt-365) & # only look back one year 
                  !is.na(date_diag)) |> 
  dplyr::group_by(lopnr, diagnosis) |>
  dplyr::arrange(date_diag) |>
  dplyr::slice(1) |>
  dplyr::select(lopnr, diagnosis) |>
  dplyr::rename(comorb_1y = diagnosis)

# cancer (excluding non-melanoma skin cancer) in previous year
ICD_cancer <- c(
  paste0("C0", 1:9),
  paste0("C", c(10:26, 30:34, 37:41, 43, 45:58, 60:76, 81:86, 88, 90:97))
)
ICD_cancer <- paste(paste0("^", ICD_cancer), collapse = "|")

# extract cancer diagnoses
diagnosis_cancer <- diagnoses_per_ID_1y |>
  dplyr::mutate(cancer = grepl(ICD_cancer, comorb_1y)) |>
  dplyr::group_by(lopnr) |>
  dplyr::summarise(cancer = as.integer(any(cancer))) # keep one row per lopnr

# add cancer indicator to cohort
cohort <- cohort |>
  dplyr::left_join(diagnosis_cancer, by="lopnr") |>
  dplyr::mutate(cancer = tidyr::replace_na(cancer, 0))

################################################################################
### Add medications of past half a year to cohort ##############################
################################################################################
# obtain index_dt from cohort
medication_index_dt <- medication |>
  dplyr::select(lopnr, atc, edatum) |>
  dplyr::left_join(cohort, by="lopnr") |>
  dplyr::select(lopnr, atc, edatum, index_dt)

# select medication per ID, ordered on date
medication_per_ID <- medication_index_dt |>
  dplyr::filter(edatum >= (index_dt-183) & # only look back half a year
                  edatum <= index_dt &     # until index date
                  !is.na(edatum)) |>
  dplyr::group_by(lopnr, atc) |>
  dplyr::arrange(edatum) |>
  dplyr::slice(1) |>
  dplyr::select(lopnr, atc) |>
  dplyr::rename(medication = atc)

# Obtain medications
medications_index <- data.frame(names=c("bblock",
                                      "hypoglycemic",
                                      "ccb",
                                      "diur",
                                      "rasi",
                                      "lipid",
                                      "nsaid"),
                              ATC=c("^C07",  # beta blockers
                                    "^A10",  # hypoglycemic agents
                                    "^C08",  # calcium channel blockers
                                    "^C03",  # diuretics
                                    "^C09A|^C09B|^C09C|^C09D", # RASi
                                    "^C10",  # statins
                                    "^M01A")) # NSAID
for (ATC_index in 1:length(medications_index$names)){
  # indicate for each ID if the ICD was present
  medication_per_ATC <- medication_per_ID |>
    dplyr::mutate(ATC_present = grepl(medications_index$ATC[ATC_index], medication)) |>
    dplyr::group_by(lopnr) |>
    dplyr::summarise(ATC_present = as.integer(any(ATC_present))) |> # keep one row per lopnr
    dplyr::rename(!!medications_index$names[ATC_index] := ATC_present)   # rename column to current ICD

  # join with cohort
  cohort <- cohort |>
    dplyr::left_join(medication_per_ATC, by="lopnr") |>
    dplyr::mutate(!!medications_index$names[ATC_index] := # those who are not in diagnoses_per_ICD did not have that comorbidity
                    tidyr::replace_na(!!dplyr::sym(medications_index$names[ATC_index]), 0))
}

# look at proportions
colMeans(cohort[17:ncol(cohort)])*100

# save updated cohort dataframe
save(cohort, file = "cohort_covariates.Rdata")
