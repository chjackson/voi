################################################################################
#### Model Inputs from Data for the Chemotherapy Model
################################################################################
## Read Data from File
dat_side_effs_SC <- read.csv("01_data_raw/Side_Effects_SC.csv",
                             header = TRUE)

## Total Number of Patients in Observed Data
n_patients <- nrow(dat_side_effs_SC) 

## Number of Patients with side effects
n_side_effects <- sum(dat_side_effs_SC$Side_Effects)

## Number of Hospitalised Patients
n_hospitalised <- sum(dat_side_effs_SC$Hosp)

## Number of Patients who Died
n_died <- sum(dat_side_effs_SC$Death)
