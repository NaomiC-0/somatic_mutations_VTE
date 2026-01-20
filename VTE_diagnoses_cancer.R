# R script for getting VTE diagnoses from the LabKey API
# Naomi Cornish
# January 2025 

# Config ####
# First run the config.R script to set up libraries and the labkey_sql function
setwd(Sys.getenv("wd"))
source('./scripts/config421.R')

# Version ####
# check which datarelease being used
print(version)
# see data release v19 notes for relevant censorship dates

# Cancer cohort #####
# load in the cancer cohort file
cancer_prog <- fread('./data/pcs_imd_stage_cancer_ncras_treatment.csv')
# list of cancer_prog participants
participant <- toString(cancer_prog$participant_id)

# ICD10 VTE codes ####
VTE_ICD10_descriptions <-fread('./data/code_lists/VTE_ICD10_algorithm.csv') 
# see supplementary tables in manuscript for details of this algorithm
VTE_ICD10 <- VTE_ICD10_descriptions %>% pull(ICD10_code)
# Codes obtained from HDRUK phenotype library PH338 and PH71 ICD10 codes (see online)
# additional codes I808-I809 and I828-829 added after protocol peer review

# HES data Notes here: https://re-docs.genomicsengland.co.uk/general_clinical/
# 'admidate' = date of each admission; 'diag_count'= number of diagnoses for that admission;
# disdate = discharge date
# 'diag_all' = concatenated (all) diagnoses for that admission

# Inpatient (hes_apc) diagnoses #####
# SQL for the INPATIENT data (hes_apc)
temp <- list()

# construct a separate SQL query for each ICD10 code. Look through the HES records for all the cancer participants
for (ICD10code in VTE_ICD10) {
sqlstr <- paste("SELECT participant_id, diag_all, admidate ",
                "FROM hes_apc",
                " WHERE participant_id IN (", participant, ")",
                " AND diag_all LIKE '%", ICD10code, "%'",
                sep = "")
temp[[ICD10code]] <- labkey_to_df(sqlstr, version, 100000) %>% as.list
 }
# this outputs a list with all the VTE diagnoses
# convert to a dataframe
list_to_df <- function(lst) {
  df <- data.frame(
    participant_id = lst[[1]],
    diag_all = lst[[2]],
    admidate = lst[[3]]
  )
  return(df)
}

all_VTE_apc <- lapply(temp, list_to_df) %>% bind_rows

# For each participant, select the first/earliest VTE diagnosis
first_VTE_temp <- all_VTE_apc %>%
  group_by(participant_id) %>%
  filter(admidate == min(admidate)) %>%
  ungroup()

# the diagnosis column is a string of multiple ICD codes (varying lengths)
#separate these into separate columns
max_splits <- max(sapply(strsplit(first_VTE_temp$diag_all, "\\|"), length)) 
# to determine how many columns required

# split the diag_all column
first_VTE_apc <- first_VTE_temp %>% 
  separate(diag_all, into = paste0("diag", 1:max_splits), sep = "\\|", fill = "right") %>%
  # note all inpatient diagnoses for each admission are listed in these tables. 
  # Select only the ICD10 codes for VTE
  # first pivot_longer to make this easier
  pivot_longer(cols = starts_with("diag"), names_to = "diag_no", values_to = "diagnosis") %>%
  # select relevant ICD10 codes
  filter(str_detect(diagnosis, paste(VTE_ICD10, collapse = "|"))) %>%
  # remove the 'diag_no' column as it doesn't matter which position the VTE diagnosis was in the list of ICD10 codes
  # all codes relate to the same 'admidate'
  select(-c(diag_no)) %>%
  #remove duplicates
  distinct

# clear environment of variables which are no longer required
rm(all_VTE_apc, max_splits)

#to clean this
first_VTE_apc_clean <- first_VTE_apc %>% 
  mutate(diagnosis = str_replace(diagnosis, "([A-Za-z]\\d{2,3}).*", "\\1"))%>%
  # check again for duplicates
  distinct

first_VTE_apc_clean2 <- first_VTE_apc_clean %>% 
  group_by(participant_id) %>%
  summarize(vte_diag = paste(unique(diagnosis), collapse = ", "),
            admidate = first(admidate))

# remove redundant environment variables ready for next workflow
rm(list=setdiff(ls(), c('first_VTE_apc_clean2', 'VTE_ICD10', 'participant', 'NONcancer', 'version', '%!in%', 'labkey_to_df')))

# Outpatient (hes_op) diagnoses #####
# SQL for the OUTPATIENT data (hes_op)
temp <- list()

# construct a separate SQL query for each ICD10 code. Look through the HES records for all the cancer participants
for (ICD10code in VTE_ICD10) {
  sqlstr <- paste("SELECT participant_id, diag_all, apptdate ",
                  "FROM hes_op",
                  " WHERE participant_id IN (", participant, ")",
                  " AND diag_all LIKE '%", ICD10code, "%'",
                  sep = "")
  temp[[ICD10code]] <- labkey_to_df(sqlstr, version, 100000) %>% as.list
  }


# this outputs a list with all the VTE diagnoses
# convert to a dataframe
list_to_df <- function(lst) {
  df <- data.frame(
    participant_id = lst[[1]],
    diag_all = lst[[2]],
    apptdate = lst[[3]]
  )
  return(df)
}

all_VTE_op <- lapply(temp, list_to_df) %>% bind_rows

# each participant appears in the df multiple times. For each participant, select the first/earliest VTE diagnosis
first_VTE_temp <- all_VTE_op %>%
  group_by(participant_id) %>%
  filter(apptdate == min(apptdate)) %>%
  ungroup()

# the diagnosis column is a string of multiple ICD codes (varying lengths)
#separate these into separate columns
max_splits <- max(sapply(strsplit(first_VTE_temp$diag_all, "\\|"), length)) 
# to determine how many columns required


# split the diag_all column
first_VTE_op <- first_VTE_temp %>% 
  separate(diag_all, into = paste0("diag", 1:(max_splits+1)), sep = "\\|", fill = "right") %>%
  # Select only the ICD10 codes for VTE
  # first pivot_longer
  pivot_longer(cols = starts_with("diag"), names_to = "diag_no", values_to = "diagnosis") %>%
  # select relevant ICD10 codes
  filter(str_detect(diagnosis, paste(VTE_ICD10, collapse = "|"))) %>%
  # all codes relate to the same 'admidate'
  select(-c(diag_no)) %>%
  #remove duplicates
  distinct

# clear environment of variables which are no longer required
rm(all_VTE_op, max_splits)

#to clean this
first_VTE_op_clean <- first_VTE_op %>% 
  mutate(diagnosis = str_replace(diagnosis, "([A-Za-z]\\d{2,3}).*", "\\1"))%>%
  # check again for duplicates
  distinct

# Note some participants are still appearing twice if they had multiple VTE ICD10 codes recorded on same date
first_VTE_op_clean %>% 
  filter(duplicated(participant_id)|duplicated(participant_id, fromLast = T))

first_VTE_op_clean2 <- first_VTE_op_clean %>% 
  group_by(participant_id) %>%
  summarize(vte_diag = paste(unique(diagnosis), collapse = ", "),
            apptdate = first(apptdate))

# remove redundant environment variables ready for next workflow
rm(list=setdiff(ls(), c('first_VTE_apc_clean2', 'first_VTE_op_clean2', 'VTE_ICD10', 'participant', 'NONcancer', 'version', '%!in%', 'labkey_to_df')))

# A+E (hes_ae) diagnoses (from 1/4/2007 to 31/3/2020) #####
# SQL for the EMERGENCY DEPT data (hes_op)
temp <- list()

# construct a separate SQL query for each ICD10 code. Look through the HES records for all the cancer participants
for (ICD10code in VTE_ICD10) {
  sqlstr <- paste("SELECT participant_id, diag_all, arrivaldate ",
                  "FROM hes_ae",
                  " WHERE participant_id IN (", participant, ")",
                  " AND diag_all LIKE '%", ICD10code, "%'",
                  sep = "")
  temp[[ICD10code]] <- labkey_to_df(sqlstr, version, 100000) %>% as.list
}


# this outputs a list with all the VTE diagnoses
# convert to a dataframe
list_to_df <- function(lst) {
  df <- data.frame(
    participant_id = lst[[1]],
    diag_all = lst[[2]],
    arrivaldate = lst[[3]]
  )
  return(df)
}

all_VTE_ae <- lapply(temp, list_to_df) %>% bind_rows

first_VTE_temp <- all_VTE_ae %>%
  group_by(participant_id) %>%
  filter(arrivaldate == min(arrivaldate)) %>%
  ungroup()

# the diagnosis column is a string of multiple ICD codes (varying lengths)
# separate these into separate columns
max_splits <- max(sapply(strsplit(first_VTE_temp$diag_all, "\\|"), length)) 
# to determine how many columns required

# split the diag_all column
first_VTE_ae <- first_VTE_temp %>% 
  separate(diag_all, into = paste0("diag", 1:(max_splits+1)), sep = "\\|", fill = "right") %>%
  # note +1 added to max_splits because I was getting  a warning message that one split was missing
  pivot_longer(cols = starts_with("diag"), names_to = "diag_no", values_to = "diagnosis") %>%
  # select relevant ICD10 codes
  filter(str_detect(diagnosis, paste(VTE_ICD10, collapse = "|"))) %>%
  # all codes relate to the same 'arrivaldate'
  select(-c(diag_no)) %>%
  #remove duplicates
  distinct

# clear environment of variables which are no longer required
rm(all_VTE_ae, max_splits)

first_VTE_ae_clean <- first_VTE_ae %>% 
  mutate(diagnosis = str_replace(diagnosis, "([A-Za-z]\\d{2,3}).*", "\\1"))%>%
  # check again for duplicates
  distinct

table(first_VTE_ae_clean$diagnosis)

first_VTE_ae_clean2 <- first_VTE_ae_clean %>% 
  group_by(participant_id) %>%
  summarize(vte_diag = paste(unique(diagnosis), collapse = ", "),
            arrivaldate = first(arrivaldate))


# remove redundant environment variables ready for next workflow
rm(list=setdiff(ls(), c('first_VTE_apc_clean2', 'first_VTE_op_clean2', 'first_VTE_ae_clean2', 'VTE_ICD10', 'participant', 'NONcancer', 'version', '%!in%', 'labkey_to_df')))

# ECDS records (A+E from 5/4/2017 to 4/8/2022) #######
# SNOMED codes not ICD10 codes
# to facilitate conversion use GEMINI code lists
# see https://github.com/GEMINI-multimorbidity

## SNOMED VTE list ####
gemini <- fread('./data/code_lists/complete_gemini_codelist.txt')
# note this code list contains codes for arterial and venous thromboembolism
vte_snomed <- gemini %>% 
  filter(condition == 'thromboembolic' &
         vocab_id == 'SNOMEDCT' &
          # remove codes with strings relating to arterial disease 
         str_detect(description, 'arter', negate =T) &
         str_detect(description, 'aorta', negate =T))

# check no superficial venous thrombosis diagnoses are included
vte_snomed %>% filter(str_detect(description, 'superficial'))

snomed_codes <- vte_snomed %>% select(code) %>% unlist

## ecds SQL query ######
temp <- list()

# construct a separate SQL query for each ICD10 code. 
for (snomed_code in snomed_codes) {
  sqlstr <- paste("SELECT participant_id, diagnosis_code_all, arrival_date ",
                  "FROM ecds",
                  " WHERE participant_id IN (", participant, ")",
                  " AND diagnosis_code_all LIKE '%", snomed_code, "%'",
                  sep = "")
  temp[[snomed_code]] <- labkey_to_df(sqlstr, version, 100000) %>% as.list
}


# this outputs a list with all the VTE diagnoses
# convert to a dataframe
list_to_df <- function(lst) {
  df <- data.frame(
    participant_id = lst[[1]],
    diag_all = lst[[2]],
    arrivaldate = lst[[3]]
  )
  return(df)
}

all_VTE_ecds <- lapply(temp, list_to_df) %>% bind_rows

# each participant appears in the df multiple times. For each participant, select the first/earliest VTE diagnosis
first_VTE_temp <- all_VTE_ecds %>%
  group_by(participant_id) %>%
  filter(arrivaldate == min(arrivaldate)) %>%
  ungroup()

# the diagnosis column is a string of multiple ICD codes (varying lengths) - separate these into separate columns
max_splits <- max(sapply(strsplit(first_VTE_temp$diag_all, "\\|"), length)) # to determine how many columns required


# split the diag_all column
first_VTE_ecds <- first_VTE_temp %>% 
  separate(diag_all, into = paste0("diag", 1:(max_splits+1)), sep = "\\|", fill = "right") %>%
  # first pivot_longer to make this easier
  pivot_longer(cols = starts_with("diag"), names_to = "diag_no", values_to = "diagnosis") %>%
  # select relevant SNOMED codes
  filter(str_detect(diagnosis, paste(snomed_codes, collapse = "|"))) %>%
  # all codes relate to the same 'arrivaldate'
  select(-c(diag_no)) %>%
  #remove duplicates
  distinct

# clear environment of variables which are no longer required
rm(all_VTE_ecds, max_splits)

first_VTE_ecds_clean <- first_VTE_ecds %>% 
  group_by(participant_id) %>%
  summarize(vte_diag = paste(unique(diagnosis), collapse = ", "),
            arrivaldate = first(arrivaldate))

# Merge HES records ####
# combine the VTE INPATIENT, OUTPATIENT, A+E diagnoses

# read in the files and rename the date column so it is consistent across all files
inpt <- first_VTE_apc_clean2 %>% rename('date_vte' = 'admidate')
outpt <- first_VTE_op_clean2 %>% rename('date_vte' = 'apptdate')
ae <- first_VTE_ae_clean2 %>% rename('date_vte' = 'arrivaldate')
ecds <- first_VTE_ecds_clean %>% rename('date_vte' = 'arrivaldate') 
rm('first_VTE_apc_clean2', 'first_VTE_op_clean2', 'first_VTE_ae_clean2', 'first_VTE_ecds_clean')

# merge the records
vte_hes <- rbind(inpt, outpt, ae, ecds) %>%
# select the first recorded date of VTE across all the records for each participant  
  group_by(participant_id) %>%
  # retain any entries from the same date
  slice_min(date_vte, with_ties = T) %>%
  # remove duplicates (so for example if I260 was recorded in hes_op and hes_apc on the same date this will be resolved by selecting unique entries)
  ungroup() %>% distinct
            
# concatenate SNOMED and ICD10 codes from same date
vte_hes_unique <- vte_hes %>% 
  group_by(participant_id) %>%
  summarize(diag_vte = paste(unique(vte_diag), collapse = ", "),
            date_vte = first(date_vte))

fwrite(vte_hes_unique, './data/hes_combined_VTE_cancercohort.csv')

# remove redundant environment variables ready for next workflow
rm(list=ls())

