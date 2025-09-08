################################################################################
## This file contains functions calculate eGFR according to different equations#
## Author: Malou Magnani #######################################################
################################################################################

################################################################################
### 1. Creatinine equations ####################################################
################################################################################
# 1.1 CKD-EPIcr2009
ckd_epi_2009_cr <- function(creatinine, age, female){
  k <- ifelse(female == 1, 62, 80) # for umol/L
  alpha <- ifelse(female == 1, -0.329, -0.411)
  return(ifelse(female == 1,
                141 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-1.209)) * (0.9929 ^ age) * 1.018,
                141 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-1.209)) * (0.9929 ^ age)))}

# 1.2 CKD-EPIcr2021
ckd_epi_2021_cr <- function(creatinine, age, female){
  k <- ifelse(female == 1, 62, 80)
  alpha <- ifelse(female == 1, -0.241, -0.302)
  return(ifelse(female == 1,
                142 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-1.200)) * (0.9938 ^ age) * 1.012,
                142 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-1.200)) * (0.9938 ^ age)))
}

################################################################################
### 2. Cystatin C equations  ####################################################
################################################################################
# 2.1 CKD-EPIcys2012
ckd_epi_2012_cys <- function(cystatin, age, female){
  return(ifelse(female == 1,
                133 * (pmin(cystatin/0.8, 1) ^ (-0.499)) * (pmax(cystatin/0.8, 1) ^ (-1.328)) * (0.9962 ^ age) * 0.932,
                133 * (pmin(cystatin/0.8, 1) ^ (-0.499)) * (pmax(cystatin/0.8, 1) ^ (-1.328)) * (0.9962 ^ age)))
}

################################################################################
### 3. Creatinine/cystatin C equations ##########################################
################################################################################
# 3.1 CKD-EPI2012crcys
ckd_epi_2012_cr_cys <- function(creatinine, cystatin, age, female){
  k <- ifelse(female == 1, 62, 80)
  alpha <- ifelse(female == 1, -0.248, -0.207)
  return(ifelse(female == 1,
                135 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-0.601)) * (pmin(cystatin/0.8, 1) ^ (-0.375)) * (pmax(cystatin/0.8, 1) ^ (-0.711)) * (0.9952 ^ age) * 0.969,
                135 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-0.601)) * (pmin(cystatin/0.8, 1) ^ (-0.375)) * (pmax(cystatin/0.8, 1) ^ (-0.711)) * (0.9952 ^ age)))
}

# 3.2 CKD-EPI2021crcys
ckd_epi_2021_cr_cys <- function(creatinine, cystatin, age, female){
  k <- ifelse(female == 1, 62, 80)
  alpha <- ifelse(female == 1, -0.219, -0.144)
  return(ifelse(female == 1,
                135 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-0.544)) * (pmin(cystatin/0.8, 1) ^ (-0.323)) * (pmax(cystatin/0.8, 1) ^ (-0.778)) * (0.9961 ^ age) * 0.963,
                135 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-0.544)) * (pmin(cystatin/0.8, 1) ^ (-0.323)) * (pmax(cystatin/0.8, 1) ^ (-0.778)) * (0.9961 ^ age)))
}