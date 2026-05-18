################################################################################
## This script is intended to load data ########################################
## Author: Malou Magnani #######################################################
################################################################################
# remove history
rm(list=ls(all.names=TRUE))

# set seed for reproducibility
set.seed(27)

# set directory to load and save data
setwd("~/Data/")

# load libraries
library(DBI)       # database connection
library(odbc)      # to load data from Kosmos
library(dplyr)     # manipulate dataframe
library(dbplyr)    # databases backend

# 0. load the necessary data tables
con <- DBI::dbConnect(odbc::odbc(), "kosmos")
tables <- dplyr::tbl(con, dbplyr::in_catalog("SCREAM2", 
                                      "INFORMATION_SCHEMA", 
                                      "TABLES")) |>
  dplyr::collect()

## 1. collect creatinine and cystatine C
catalog <- "SCREAM2"
schema <- "SCREAM2_CLEAN"
lab_values <- dplyr::tbl(con, dbplyr::in_catalog(catalog, schema, "lab_values")) |> 
  dplyr::filter(test == "crea" | test == "cystC") |>
  dplyr::collect()
save(lab_values, file = "lab_values.Rda")

## 2. collect demographics
demo <- dplyr::tbl(con, dbplyr::in_catalog(catalog, schema, "scream2_demo")) |> 
  dplyr::collect()
save(demo, file = "demographics.Rda")

## 3. collect renal replacement therapy (RRT)
rrt <- dplyr::tbl(con, dbplyr::in_catalog(catalog, schema, "snr_rrt")) |> 
  dplyr::collect()
save(rrt, file = "rrt.Rda")

## 4. collect information on mortality
death <- dplyr::tbl(con, dbplyr::in_catalog(catalog, schema, "scream2_death")) |> 
  dplyr::collect()
save(death, file = "death.Rda")

## 5. collect information on emigration
emigration <- dplyr::tbl(con, dbplyr::in_catalog(catalog, schema, "scream2_flytt")) |> 
  dplyr::collect()
save(emigration, file = "emigration.Rda")

## 6.1 collect diagnoses in hospital care
diagnoses_slv <- dplyr::tbl(con, in_catalog(catalog, schema, "slv_diagnoses_long")) |> 
  filter(diagnosis %like% "I200%" |
  diagnosis %like% "I21%" |
  diagnosis %like% "I22%" |
  diagnosis %like% "I201%" |
  diagnosis %like% "I208%" |
  diagnosis %like% "I209%" |
  diagnosis %like% "I24%" |
  diagnosis %like% "I25%" |
  diagnosis %like% "I10%" |
  diagnosis %like% "I11%" |
  diagnosis %like% "I12%" |
  diagnosis %like% "I13%" |
  diagnosis %like% "I14%" |
  diagnosis %like% "I15%" |
  diagnosis %like% "I110%" |
  diagnosis %like% "I130%" |
  diagnosis %like% "I132%" |
  diagnosis %like% "I50%" |
  diagnosis %like% "I60%" |
  diagnosis %like% "I61%" |
  diagnosis %like% "I62%" |
  diagnosis %like% "I63%" |
  diagnosis %like% "I64%" |
  diagnosis %like% "I693%" |
  diagnosis %like% "I698%" |
  diagnosis %like% "I694%" |
  diagnosis %like% "I65%" |
  diagnosis %like% "I66%" |
  diagnosis %like% "I67%" |
  diagnosis %like% "I68%" |
  diagnosis %like% "I69%" |
  diagnosis %like% "G450%" |
  diagnosis %like% "G451%" |
  diagnosis %like% "G452%" |
  diagnosis %like% "G453%" |
  diagnosis %like% "G458%" |
  diagnosis %like% "G459%" |
  diagnosis %like% "G46%" |
  diagnosis %like% "I44%" |
  diagnosis %like% "I45%" |
  diagnosis %like% "I46%" |
  diagnosis %like% "I47%" |
  diagnosis %like% "I48%" |
  diagnosis %like% "I49%" |
  diagnosis %like% "I70%" |
  diagnosis %like% "I72%" |
  diagnosis %like% "I73%" |
  diagnosis %like% "E10%" |
  diagnosis %like% "E11%" |
  diagnosis %like% "E12%" |
  diagnosis %like% "E13%" |
  diagnosis %like% "E14%" |
  diagnosis %like% "J44%" |
  diagnosis %like% "C01%" |
  diagnosis %like% "C02%" |
  diagnosis %like% "C03%" |
  diagnosis %like% "C04%" |
  diagnosis %like% "C05%" |
  diagnosis %like% "C06%" |
  diagnosis %like% "C07%" |
  diagnosis %like% "C08%" |
  diagnosis %like% "C09%" |
  diagnosis %like% "C10%" |
  diagnosis %like% "C11%" |
  diagnosis %like% "C12%" |
  diagnosis %like% "C13%" |
  diagnosis %like% "C14%" |
  diagnosis %like% "C15%" |
  diagnosis %like% "C16%" |
  diagnosis %like% "C17%" |
  diagnosis %like% "C18%" |
  diagnosis %like% "C19%" |
  diagnosis %like% "C20%" |
  diagnosis %like% "C21%" |
  diagnosis %like% "C22%" |
  diagnosis %like% "C23%" |
  diagnosis %like% "C24%" |
  diagnosis %like% "C25%" |
  diagnosis %like% "C26%" |
  diagnosis %like% "C30%" |
  diagnosis %like% "C31%" |
  diagnosis %like% "C32%" |
  diagnosis %like% "C33%" |
  diagnosis %like% "C34%" |
  diagnosis %like% "C37%" |
  diagnosis %like% "C38%" |
  diagnosis %like% "C39%" |
  diagnosis %like% "C40%" |
  diagnosis %like% "C41%" |
  diagnosis %like% "C43%" |
  diagnosis %like% "C45%" |
  diagnosis %like% "C46%" |
  diagnosis %like% "C47%" |
  diagnosis %like% "C48%" |
  diagnosis %like% "C49%" |
  diagnosis %like% "C50%" |
  diagnosis %like% "C51%" |
  diagnosis %like% "C52%" |
  diagnosis %like% "C53%" |
  diagnosis %like% "C54%" |
  diagnosis %like% "C55%" |
  diagnosis %like% "C56%" |
  diagnosis %like% "C57%" |
  diagnosis %like% "C58%" |
  diagnosis %like% "C60%" |
  diagnosis %like% "C61%" |
  diagnosis %like% "C62%" |
  diagnosis %like% "C63%" |
  diagnosis %like% "C64%" |
  diagnosis %like% "C65%" |
  diagnosis %like% "C66%" |
  diagnosis %like% "C67%" |
  diagnosis %like% "C68%" |
  diagnosis %like% "C69%" |
  diagnosis %like% "C70%" |
  diagnosis %like% "C71%" |
  diagnosis %like% "C72%" |
  diagnosis %like% "C73%" |
  diagnosis %like% "C74%" |
  diagnosis %like% "C75%" |
  diagnosis %like% "C76%" |
  diagnosis %like% "C81%" |
  diagnosis %like% "C82%" |
  diagnosis %like% "C83%" |
  diagnosis %like% "C84%" |
  diagnosis %like% "C85%" |
  diagnosis %like% "C86%" |
  diagnosis %like% "C88%" |
  diagnosis %like% "C90%" |
  diagnosis %like% "C91%" |
  diagnosis %like% "C92%" |
  diagnosis %like% "C93%" |
  diagnosis %like% "C94%" |
  diagnosis %like% "C95%" |
  diagnosis %like% "C96%" |
  diagnosis %like% "C97%" |
  diagnosis %like% "B18%" |
  diagnosis %like% "I850%" |
  diagnosis %like% "I859%" |
  diagnosis %like% "I982%" |
  diagnosis %like% "K70%" |
  diagnosis %like% "K71%" |
  diagnosis %like% "K72%" |
  diagnosis %like% "K73%" |
  diagnosis %like% "K74%" |
  diagnosis %like% "K75%" |
  diagnosis %like% "K76%" |
  diagnosis %like% "K77%" |
  diagnosis %like% "Z524%") |> 
  collect()
save(diagnoses_slv, file = "diagnoses_slv.Rda")

## 6.2 collect diagnoses in outpatient care
diagnoses_ovr <- tbl(con, in_catalog(catalog, schema, "ovr_diagnoses_long")) |> 
  filter(diagnosis %like% "I200%" |
           diagnosis %like% "I21%" |
           diagnosis %like% "I22%" |
           diagnosis %like% "I201%" |
           diagnosis %like% "I208%" |
           diagnosis %like% "I209%" |
           diagnosis %like% "I24%" |
           diagnosis %like% "I25%" |
           diagnosis %like% "I10%" |
           diagnosis %like% "I11%" |
           diagnosis %like% "I12%" |
           diagnosis %like% "I13%" |
           diagnosis %like% "I14%" |
           diagnosis %like% "I15%" |
           diagnosis %like% "I110%" |
           diagnosis %like% "I130%" |
           diagnosis %like% "I132%" |
           diagnosis %like% "I50%" |
           diagnosis %like% "I60%" |
           diagnosis %like% "I61%" |
           diagnosis %like% "I62%" |
           diagnosis %like% "I63%" |
           diagnosis %like% "I64%" |
           diagnosis %like% "I693%" |
           diagnosis %like% "I698%" |
           diagnosis %like% "I694%" |
           diagnosis %like% "I65%" |
           diagnosis %like% "I66%" |
           diagnosis %like% "I67%" |
           diagnosis %like% "I68%" |
           diagnosis %like% "I69%" |
           diagnosis %like% "G450%" |
           diagnosis %like% "G451%" |
           diagnosis %like% "G452%" |
           diagnosis %like% "G453%" |
           diagnosis %like% "G458%" |
           diagnosis %like% "G459%" |
           diagnosis %like% "G46%" |
           diagnosis %like% "I44%" |
           diagnosis %like% "I45%" |
           diagnosis %like% "I46%" |
           diagnosis %like% "I47%" |
           diagnosis %like% "I48%" |
           diagnosis %like% "I49%" |
           diagnosis %like% "I70%" |
           diagnosis %like% "I72%" |
           diagnosis %like% "I73%" |
           diagnosis %like% "E10%" |
           diagnosis %like% "E11%" |
           diagnosis %like% "E12%" |
           diagnosis %like% "E13%" |
           diagnosis %like% "E14%" |
           diagnosis %like% "J44%" |
           diagnosis %like% "C01%" |
           diagnosis %like% "C02%" |
           diagnosis %like% "C03%" |
           diagnosis %like% "C04%" |
           diagnosis %like% "C05%" |
           diagnosis %like% "C06%" |
           diagnosis %like% "C07%" |
           diagnosis %like% "C08%" |
           diagnosis %like% "C09%" |
           diagnosis %like% "C10%" |
           diagnosis %like% "C11%" |
           diagnosis %like% "C12%" |
           diagnosis %like% "C13%" |
           diagnosis %like% "C14%" |
           diagnosis %like% "C15%" |
           diagnosis %like% "C16%" |
           diagnosis %like% "C17%" |
           diagnosis %like% "C18%" |
           diagnosis %like% "C19%" |
           diagnosis %like% "C20%" |
           diagnosis %like% "C21%" |
           diagnosis %like% "C22%" |
           diagnosis %like% "C23%" |
           diagnosis %like% "C24%" |
           diagnosis %like% "C25%" |
           diagnosis %like% "C26%" |
           diagnosis %like% "C30%" |
           diagnosis %like% "C31%" |
           diagnosis %like% "C32%" |
           diagnosis %like% "C33%" |
           diagnosis %like% "C34%" |
           diagnosis %like% "C37%" |
           diagnosis %like% "C38%" |
           diagnosis %like% "C39%" |
           diagnosis %like% "C40%" |
           diagnosis %like% "C41%" |
           diagnosis %like% "C43%" |
           diagnosis %like% "C45%" |
           diagnosis %like% "C46%" |
           diagnosis %like% "C47%" |
           diagnosis %like% "C48%" |
           diagnosis %like% "C49%" |
           diagnosis %like% "C50%" |
           diagnosis %like% "C51%" |
           diagnosis %like% "C52%" |
           diagnosis %like% "C53%" |
           diagnosis %like% "C54%" |
           diagnosis %like% "C55%" |
           diagnosis %like% "C56%" |
           diagnosis %like% "C57%" |
           diagnosis %like% "C58%" |
           diagnosis %like% "C60%" |
           diagnosis %like% "C61%" |
           diagnosis %like% "C62%" |
           diagnosis %like% "C63%" |
           diagnosis %like% "C64%" |
           diagnosis %like% "C65%" |
           diagnosis %like% "C66%" |
           diagnosis %like% "C67%" |
           diagnosis %like% "C68%" |
           diagnosis %like% "C69%" |
           diagnosis %like% "C70%" |
           diagnosis %like% "C71%" |
           diagnosis %like% "C72%" |
           diagnosis %like% "C73%" |
           diagnosis %like% "C74%" |
           diagnosis %like% "C75%" |
           diagnosis %like% "C76%" |
           diagnosis %like% "C81%" |
           diagnosis %like% "C82%" |
           diagnosis %like% "C83%" |
           diagnosis %like% "C84%" |
           diagnosis %like% "C85%" |
           diagnosis %like% "C86%" |
           diagnosis %like% "C88%" |
           diagnosis %like% "C90%" |
           diagnosis %like% "C91%" |
           diagnosis %like% "C92%" |
           diagnosis %like% "C93%" |
           diagnosis %like% "C94%" |
           diagnosis %like% "C95%" |
           diagnosis %like% "C96%" |
           diagnosis %like% "C97%" |
           diagnosis %like% "B18%" |
           diagnosis %like% "I850%" |
           diagnosis %like% "I859%" |
           diagnosis %like% "I982%" |
           diagnosis %like% "K70%" |
           diagnosis %like% "K71%" |
           diagnosis %like% "K72%" |
           diagnosis %like% "K73%" |
           diagnosis %like% "K74%" |
           diagnosis %like% "K75%" |
           diagnosis %like% "K76%" |
           diagnosis %like% "K77%" |
           diagnosis %like% "Z524%") |> 
  dplyr::collect()
save(diagnoses_ovr, file = "diagnoses_ovr.Rda")

## 6.3 collect diagnoses in primary care visit
diagnoses_kon <- tbl(con, in_catalog(catalog, schema, "kon_diagnoses_long")) |> 
  filter(diagnosis %like% "I200%" |
           diagnosis %like% "I21%" |
           diagnosis %like% "I22%" |
           diagnosis %like% "I201%" |
           diagnosis %like% "I208%" |
           diagnosis %like% "I209%" |
           diagnosis %like% "I24%" |
           diagnosis %like% "I25%" |
           diagnosis %like% "I10%" |
           diagnosis %like% "I11%" |
           diagnosis %like% "I12%" |
           diagnosis %like% "I13%" |
           diagnosis %like% "I14%" |
           diagnosis %like% "I15%" |
           diagnosis %like% "I110%" |
           diagnosis %like% "I130%" |
           diagnosis %like% "I132%" |
           diagnosis %like% "I50%" |
           diagnosis %like% "I60%" |
           diagnosis %like% "I61%" |
           diagnosis %like% "I62%" |
           diagnosis %like% "I63%" |
           diagnosis %like% "I64%" |
           diagnosis %like% "I693%" |
           diagnosis %like% "I698%" |
           diagnosis %like% "I694%" |
           diagnosis %like% "I65%" |
           diagnosis %like% "I66%" |
           diagnosis %like% "I67%" |
           diagnosis %like% "I68%" |
           diagnosis %like% "I69%" |
           diagnosis %like% "G450%" |
           diagnosis %like% "G451%" |
           diagnosis %like% "G452%" |
           diagnosis %like% "G453%" |
           diagnosis %like% "G458%" |
           diagnosis %like% "G459%" |
           diagnosis %like% "G46%" |
           diagnosis %like% "I44%" |
           diagnosis %like% "I45%" |
           diagnosis %like% "I46%" |
           diagnosis %like% "I47%" |
           diagnosis %like% "I48%" |
           diagnosis %like% "I49%" |
           diagnosis %like% "I70%" |
           diagnosis %like% "I72%" |
           diagnosis %like% "I73%" |
           diagnosis %like% "E10%" |
           diagnosis %like% "E11%" |
           diagnosis %like% "E12%" |
           diagnosis %like% "E13%" |
           diagnosis %like% "E14%" |
           diagnosis %like% "J44%" |
           diagnosis %like% "C01%" |
           diagnosis %like% "C02%" |
           diagnosis %like% "C03%" |
           diagnosis %like% "C04%" |
           diagnosis %like% "C05%" |
           diagnosis %like% "C06%" |
           diagnosis %like% "C07%" |
           diagnosis %like% "C08%" |
           diagnosis %like% "C09%" |
           diagnosis %like% "C10%" |
           diagnosis %like% "C11%" |
           diagnosis %like% "C12%" |
           diagnosis %like% "C13%" |
           diagnosis %like% "C14%" |
           diagnosis %like% "C15%" |
           diagnosis %like% "C16%" |
           diagnosis %like% "C17%" |
           diagnosis %like% "C18%" |
           diagnosis %like% "C19%" |
           diagnosis %like% "C20%" |
           diagnosis %like% "C21%" |
           diagnosis %like% "C22%" |
           diagnosis %like% "C23%" |
           diagnosis %like% "C24%" |
           diagnosis %like% "C25%" |
           diagnosis %like% "C26%" |
           diagnosis %like% "C30%" |
           diagnosis %like% "C31%" |
           diagnosis %like% "C32%" |
           diagnosis %like% "C33%" |
           diagnosis %like% "C34%" |
           diagnosis %like% "C37%" |
           diagnosis %like% "C38%" |
           diagnosis %like% "C39%" |
           diagnosis %like% "C40%" |
           diagnosis %like% "C41%" |
           diagnosis %like% "C43%" |
           diagnosis %like% "C45%" |
           diagnosis %like% "C46%" |
           diagnosis %like% "C47%" |
           diagnosis %like% "C48%" |
           diagnosis %like% "C49%" |
           diagnosis %like% "C50%" |
           diagnosis %like% "C51%" |
           diagnosis %like% "C52%" |
           diagnosis %like% "C53%" |
           diagnosis %like% "C54%" |
           diagnosis %like% "C55%" |
           diagnosis %like% "C56%" |
           diagnosis %like% "C57%" |
           diagnosis %like% "C58%" |
           diagnosis %like% "C60%" |
           diagnosis %like% "C61%" |
           diagnosis %like% "C62%" |
           diagnosis %like% "C63%" |
           diagnosis %like% "C64%" |
           diagnosis %like% "C65%" |
           diagnosis %like% "C66%" |
           diagnosis %like% "C67%" |
           diagnosis %like% "C68%" |
           diagnosis %like% "C69%" |
           diagnosis %like% "C70%" |
           diagnosis %like% "C71%" |
           diagnosis %like% "C72%" |
           diagnosis %like% "C73%" |
           diagnosis %like% "C74%" |
           diagnosis %like% "C75%" |
           diagnosis %like% "C76%" |
           diagnosis %like% "C81%" |
           diagnosis %like% "C82%" |
           diagnosis %like% "C83%" |
           diagnosis %like% "C84%" |
           diagnosis %like% "C85%" |
           diagnosis %like% "C86%" |
           diagnosis %like% "C88%" |
           diagnosis %like% "C90%" |
           diagnosis %like% "C91%" |
           diagnosis %like% "C92%" |
           diagnosis %like% "C93%" |
           diagnosis %like% "C94%" |
           diagnosis %like% "C95%" |
           diagnosis %like% "C96%" |
           diagnosis %like% "C97%" |
           diagnosis %like% "B18%" |
           diagnosis %like% "I850%" |
           diagnosis %like% "I859%" |
           diagnosis %like% "I982%" |
           diagnosis %like% "K70%" |
           diagnosis %like% "K71%" |
           diagnosis %like% "K72%" |
           diagnosis %like% "K73%" |
           diagnosis %like% "K74%" |
           diagnosis %like% "K75%" |
           diagnosis %like% "K76%" |
           diagnosis %like% "K77%" |
           diagnosis %like% "Z524%") |> 
  dplyr::collect()
save(diagnoses_kon, file = "diagnoses_kon.Rda")

## 7. collect information on medication
medication <- tbl(con, in_catalog(catalog, schema, "v_x_lmed_all")) |> 
  filter(atc %like% "C07%" |
           atc %like% "A10%" |
           atc %like% "C08%" |
           atc %like% "C07%" |
           atc %like% "C03%" |
           atc %like% "C09%" |
           atc %like% "C10%" |
           atc %like% "M01A%") |> 
  dplyr::collect()
save(medication, file = "medication.Rda")