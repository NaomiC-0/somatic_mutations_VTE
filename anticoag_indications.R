# Script for identifying participants with a prior indication for anticoagulaiton or antiplatelet therapy 
# in hes_apc or hes_op and hes_ae
# N Cornish

# First run the config.R script to set up libraries and the labkey_sql function
setwd(Sys.getenv("wd"))
source('./scripts/config421.R')

# ICD10 codes of interest ####
# Manually  mapped from the list of indication for anticoag/antiplatelets in the NICE treatment summaries to ICD10 directly 
# see supplementary tables in manuscript for this list
anticoag_icd10 <- fread('./data/code_lists/anticoagulation_indications_ICD10.csv')
print(anticoag_icd10)

# create string of the relevant codes in order to pull codes from the hes tables on labkey
anticoag_indications <- anticoag_icd10 %>% select(ICD10_code) %>% unlist

# Cancer cohort #####
# load in the cohort file produced in the format_cancer_cohort script
cancer_prog <- fread('./data/pcs_imd_stage_cancer_ncras_treatment.csv')

# list of cancer_prog participants
participant <- toString(cancer_prog$participant_id)

# Notes for HES tables here: https://re-docs.genomicsengland.co.uk/general_clinical/
# 'admidate' = date of each admission; 'diag_count'= number of diagnoses for that admission;
# 'diag_all' = concatenated (all) diagnoses for that admission

# Inpatient (hes_apc) diagnoses #####
temp <- list()
# construct a separate SQL query for each ICD10 code. Look through the HES records for all the cancer participants
for (ICD10code in anticoag_indications) {
  sqlstr <- paste("SELECT participant_id, diag_all, admidate ",
                  "FROM hes_apc",
                  " WHERE participant_id IN (", participant, ")",
                  " AND diag_all LIKE '%", ICD10code, "%'",
                  sep = "")
  temp[[ICD10code]] <- labkey_to_df(sqlstr, version, 100000) %>% as.list
}
# this outputs a list with all the anticoag indications
# convert to a dataframe
list_to_df <- function(lst) {
  df <- data.frame(
    participant_id = lst[[1]],
    diag_all = lst[[2]],
    admidate = lst[[3]]
  )
  return(df)
}

all_anticoag_apc <- lapply(temp, list_to_df) %>% bind_rows

# For each participant, select the first/earliest anticoag indication
first_anticoag_temp <- all_anticoag_apc %>%
  group_by(participant_id) %>%
  filter(admidate == min(admidate)) %>%
  ungroup()

# the diagnosis column is a string of multiple ICD codes (varying lengths) - separate these into separate columns
max_splits <- max(sapply(strsplit(first_anticoag_temp$diag_all, "\\|"), length)) # to determine how many columns required

# split the diag_all column
first_anticoag_apc <- first_anticoag_temp %>% 
  separate(diag_all, into = paste0("diag", 1:max_splits), sep = "\\|", fill = "right") %>%
  # note all inpatient diagnoses for each admission are listed in these tables. 
  # Select only the ICD10 codes for anticoag
  # first pivot_longer to make this easier
  pivot_longer(cols = starts_with("diag"), names_to = "diag_no", values_to = "diagnosis") %>%
  # select relevant ICD10 codes
  filter(str_detect(diagnosis, paste(anticoag_indications, collapse = "|"))) %>%
  # remove the 'diag_no' column as it doesn't matter which position the anticoag indication was in the list of ICD10 codes
  # all codes relate to the same 'admidate'
  select(-c(diag_no)) %>%
  #remove duplicates
  distinct

# clear environment of variables which are no longer required
rm(all_anticoag_apc, max_splits)

table(first_anticoag_apc$diagnosis)

#to clean this (removing unnecessary characters such as '-' or 'X'
first_anticoag_apc_clean <- first_anticoag_apc %>% 
  mutate(diagnosis = str_replace(diagnosis, "([A-Za-z]\\d{2,3}).*", "\\1"))%>%
  # check again for duplicates
  distinct

table(first_anticoag_apc_clean$diagnosis)

# Concatenate different diagnoses relating to the same date into a single column

first_anticoag_apc_clean2 <- first_anticoag_apc_clean %>% 
  group_by(participant_id) %>%
  summarize(anticoag_diag = paste(unique(diagnosis), collapse = ", "),
            admidate = first(admidate))

first_anticoag_apc_clean2 %>% select(participant_id) %>% anyDuplicated()
# remove redundant environment variables
rm(list=setdiff(ls(), c('first_anticoag_apc_clean2',
                        'anticoag_indications', 'participant', 'NONcancer', 'version', '%!in%', 'labkey_to_df')))

# Outpatient (hes_op) diagnoses #####
# SQL for the OUTPATIENT data (hes_op)
temp <- list()

# construct a separate SQL query for each ICD10 code. Look through the HES records for all the cancer participants
for (ICD10code in anticoag_indications) {
  sqlstr <- paste("SELECT participant_id, diag_all, apptdate ",
                  "FROM hes_op",
                  " WHERE participant_id IN (", participant, ")",
                  " AND diag_all LIKE '%", ICD10code, "%'",
                  sep = "")
  temp[[ICD10code]] <- labkey_to_df(sqlstr, version, 100000) %>% as.list
}


# this outputs a list with all the anticoag indications
# convert to a dataframe
list_to_df <- function(lst) {
  df <- data.frame(
    participant_id = lst[[1]],
    diag_all = lst[[2]],
    apptdate = lst[[3]]
  )
  return(df)
}

all_anticoag_op <- lapply(temp, list_to_df) %>% bind_rows

# For each participant, select the first/earliest anticoag indication
first_anticoag_temp <- all_anticoag_op %>%
  group_by(participant_id) %>%
  filter(apptdate == min(apptdate)) %>%
  ungroup()

# the diagnosis column is a string of multiple ICD codes (varying lengths) - separate these into separate columns
max_splits <- max(sapply(strsplit(first_anticoag_temp$diag_all, "\\|"), length)) # to determine how many columns required

# split the diag_all column
first_anticoag_op <- first_anticoag_temp %>% 
  separate(diag_all, into = paste0("diag", 1:(max_splits+1)), sep = "\\|", fill = "right") %>%
  # Select only the ICD10 codes for anticoag
  # first pivot_longer to make this easier
  pivot_longer(cols = starts_with("diag"), names_to = "diag_no", values_to = "diagnosis") %>%
  # select relevant ICD10 codes
  filter(str_detect(diagnosis, paste(anticoag_indications, collapse = "|"))) %>%
  # remove the 'diag_no' column as it doesn't matter which position the anticoag indication was in the list of ICD10 codes
  # all codes relate to the same 'admidate'
  select(-c(diag_no)) %>%
  #remove duplicates
  distinct

# clear environment of variables which are no longer required
rm(all_anticoag_op, max_splits)

table(first_anticoag_op$diagnosis)

#to clean this
first_anticoag_op_clean <- first_anticoag_op %>% 
  mutate(diagnosis = str_replace(diagnosis, "([A-Za-z]\\d{2,3}).*", "\\1"))%>%
  # check again for duplicates
  distinct

table(first_anticoag_op_clean$diagnosis)

# Since the specific anticoag indication is not required: concatenate these into a single column
first_anticoag_op_clean2 <- first_anticoag_op_clean %>% 
  group_by(participant_id) %>%
  summarize(anticoag_diag = paste(unique(diagnosis), collapse = ", "),
            apptdate = first(apptdate))

# check for duplicates
first_anticoag_op_clean2 %>% select(participant_id) %>% anyDuplicated()
## save outpatient diagnoses ####

# remove redundant environment variables
rm(list=setdiff(ls(), c('first_anticoag_apc_clean2', 'first_anticoag_op_clean2',
                        'anticoag_indications', 'participant', 'NONcancer', 'version', '%!in%', 'labkey_to_df')))

# A+E (hes_ae) diagnoses  #####
# SQL for the EMERGENCY DEPT data (hes_op)
temp <- list()

# construct a separate SQL query for each ICD10 code. Look through the HES records for all the cancer participants
for (ICD10code in anticoag_indications) {
  sqlstr <- paste("SELECT participant_id, diag_all, arrivaldate ",
                  "FROM hes_ae",
                  " WHERE participant_id IN (", participant, ")",
                  " AND diag_all LIKE '%", ICD10code, "%'",
                  sep = "")
  temp[[ICD10code]] <- labkey_to_df(sqlstr, version, 100000) %>% as.list
}


# this outputs a list with all the anticoag indications
# convert to a dataframe
list_to_df <- function(lst) {
  df <- data.frame(
    participant_id = lst[[1]],
    diag_all = lst[[2]],
    arrivaldate = lst[[3]]
  )
  return(df)
}

all_anticoag_ae <- lapply(temp, list_to_df) %>% bind_rows

# each participant appears in the df multiple times. For each participant, select the first/earliest anticoag indication
first_anticoag_temp <- all_anticoag_ae %>%
  group_by(participant_id) %>%
  filter(arrivaldate == min(arrivaldate)) %>%
  ungroup()

# the diagnosis column is a string of multiple ICD codes (varying lengths) - separate these into separate columns
max_splits <- max(sapply(strsplit(first_anticoag_temp$diag_all, "\\|"), length)) # to determine how many columns required

# split the diag_all column
first_anticoag_ae <- first_anticoag_temp %>% 
  separate(diag_all, into = paste0("diag", 1:(max_splits+1)), sep = "\\|", fill = "right") %>%
  pivot_longer(cols = starts_with("diag"), names_to = "diag_no", values_to = "diagnosis") %>%
  # select relevant ICD10 codes
  filter(str_detect(diagnosis, paste(anticoag_indications, collapse = "|"))) %>%
  # remove the 'diag_no' column as it doesn't matter which position the anticoag indication was in the list of ICD10 codes
  # all codes relate to the same 'arrivaldate'
  select(-c(diag_no)) %>%
  #remove duplicates
  distinct

# clear environment of variables which are no longer required
rm(all_anticoag_ae, max_splits)

first_anticoag_ae_clean <- first_anticoag_ae %>% 
  mutate(diagnosis = str_replace(diagnosis, "([A-Za-z]\\d{2,3}).*", "\\1"))%>%
  # check again for duplicates
  distinct

table(first_anticoag_ae_clean$diagnosis)


# For now the specific anticoag indication is not required: concatenate these into a single column if participants appear more than once in the table
first_anticoag_ae_clean2 <- first_anticoag_ae_clean %>% 
  group_by(participant_id) %>%
  summarize(anticoag_diag = paste(unique(diagnosis), collapse = ", "),
            arrivaldate = first(arrivaldate))

first_anticoag_ae_clean2 %>% select(participant_id) %>% anyDuplicated()

# remove redundant environment variables ready for next workflow
rm(list=setdiff(ls(), c('first_anticoag_apc_clean2', 'first_anticoag_op_clean2', 'first_anticoag_ae_clean2',
                        'anticoag_indications', 'participant', 'NONcancer', 'version', '%!in%', 'labkey_to_df')))

 
# ECDS records (A+E from 5/4/2017 to 4/8/2022) #######
# checked the datarelease notes for version 19: censor date for ecds records is still 2022
# these use SNOMED codes not ICD10 codes
# to facilitate conversion I have downloaded the GEMINI code lists
# see https://github.com/GEMINI-multimorbidity for the citations

## SNOMED anticoag list ####
gemini <- fread('./data/code_lists/complete_gemini_codelist.txt')
# note this code list contains codes for arterial and venous thromboembolism
anticoag_snomed <- gemini %>% 
  filter(condition == 'AF'|condition == 'coronary_heart'|
           condition == 'isch_stroke'|condition == 'isch_stroke_tia'|
           condition == 'stroke_all') %>%
  filter(vocab_id == 'SNOMEDCT' &
           # remove codes with strings relating to cerebral haemorrhage
           str_detect(description, 'haemorrhage', negate =T))

# save these for use in appendix
fwrite(anticoag_snomed, './data/code_lists/anticoagulation_indications_SNOMEDcodes.csv', row.names = F)

# check only cerebral haemorrhage disease strings removed with the above code
gemini %>% 
  filter(condition == 'AF'|condition == 'coronary_heart'|
           condition == 'isch_stroke'|condition == 'isch_stroke_tia'|
           condition == 'stroke_all') %>%
  filter(vocab_id == 'SNOMEDCT' &
           # remove codes with strings relating to cerebral haemorrhage
           str_detect(description, 'haemorrhage', negate =F)) %>% View

snomed_codes <- anticoag_snomed %>% select(code) %>% unlist

## ecds SQL query ######
temp <- list()

# construct a separate SQL query for each ICD10 code. Look through the HES records for all the cancer participants
for (snomed_code in snomed_codes) {
  sqlstr <- paste("SELECT participant_id, diagnosis_code_all, arrival_date ",
                  "FROM ecds",
                  " WHERE participant_id IN (", participant, ")",
                  " AND diagnosis_code_all LIKE '%", snomed_code, "%'",
                  sep = "")
  temp[[snomed_code]] <- labkey_to_df(sqlstr, version, 100000) %>% as.list
}


# this outputs a list with all the anticoag indications
# convert to a dataframe
list_to_df <- function(lst) {
  df <- data.frame(
    participant_id = lst[[1]],
    diag_all = lst[[2]],
    arrivaldate = lst[[3]]
  )
  return(df)
}

all_anticoag_ecds <- lapply(temp, list_to_df) %>% bind_rows

first_anticoag_temp <- all_anticoag_ecds %>%
  group_by(participant_id) %>%
  filter(arrivaldate == min(arrivaldate)) %>%
  ungroup()

# the diagnosis column is a string of multiple ICD codes (varying lengths) - separate these into separate columns
max_splits <- max(sapply(strsplit(first_anticoag_temp$diag_all, "\\|"), length)) # to determine how many columns required

# split the diag_all column
first_anticoag_ecds <- first_anticoag_temp %>% 
# split the diag_all columns
  separate(diag_all, into = paste0("diag", 1:(max_splits+1)), sep = "\\|", fill = "right") %>%
  # Select only the SNOMED codes for anticoag
  # first pivot_longer to make this easier
  pivot_longer(cols = starts_with("diag"), names_to = "diag_no", values_to = "diagnosis") %>%
  # select relevant SNOMED codes
  filter(str_detect(diagnosis, paste(snomed_codes, collapse = "|"))) %>%
  # remove the 'diag_no' column as it doesn't matter which position the anticoag indication was in the list of anticoag_snomed codes
  # all codes relate to the same 'arrivaldate'
  select(-c(diag_no)) %>%
  #remove duplicates
  distinct

# clear environment of variables which are no longer required
rm(all_anticoag_ecds, max_splits)

ecds_codes <- first_anticoag_ecds %>% select(diagnosis) %>% distinct %>% unlist
# to view what these codes relate to:
gemini %>% filter(code %in% ecds_codes)
# this confirms that the codes which have been pulled from ecds data are all consistent with indications for anticoagulation


# For now the specific anticoag indication is not relevant: concatenate these into a single column
first_anticoag_ecds_clean <- first_anticoag_ecds %>% 
  group_by(participant_id) %>%
  summarize(anticoag_diag = paste(unique(diagnosis), collapse = ", "),
            arrivaldate = first(arrivaldate))

# remove redundant environment variables ready for next workflow
rm(list=setdiff(ls(), c('first_anticoag_apc_clean2', 'first_anticoag_op_clean2', 'first_anticoag_ae_clean2',
                        'first_anticoag_ecds_clean',
                        'anticoag_indications', 'participant', 'NONcancer', 'version', '%!in%', 'labkey_to_df')))

# Merge HES records ####
# combine the anticoag INPATIENT, OUTPATIENT, A+E diagnoses

# read in the files and rename the date column so it is consistent across all files
inpt <- first_anticoag_apc_clean2 %>% rename('date_anticoag' = 'admidate')
outpt <- first_anticoag_op_clean2 %>% rename('date_anticoag' = 'apptdate')
ae <- first_anticoag_ae_clean2 %>% rename('date_anticoag' = 'arrivaldate')
ecds <- first_anticoag_ecds_clean %>% rename('date_anticoag' = 'arrivaldate') 
rm(first_anticoag_apc_clean2, first_anticoag_op_clean2, first_anticoag_ae_clean2, first_anticoag_ecds_clean)

# merge the records
anticoag_hes <- rbind(inpt, outpt, ae, ecds) %>%
  # select the first recorded date of an anticoag indication across all the records for each participant  
  group_by(participant_id) %>%
  # retain any entries from the same date
  slice_min(date_anticoag, with_ties = T) %>%
  # remove duplicates (so for example if I260 was recorded in hes_op and hes_apc on the same date this will be resolved by selecting unique entries)
  ungroup() %>% distinct

# If duplicate entries from the same date persist (across the different electronic health record sources),
# concatenate the diagnoses for a single participant into one row so that each participant appears only once
anticoag_hes_unique <- anticoag_hes %>% 
  group_by(participant_id) %>%
  summarize(diag_anticoag = paste(unique(anticoag_diag), collapse = ", "),
            date_anticoag = first(date_anticoag))

# confirm no duplicates
anticoag_hes_unique %>% select(participant_id) %>% anyDuplicated()

rm(anticoag_hes)

# save
fwrite(anticoag_hes_unique, './data/hes_combined_anticoag_cancercohort.csv')


# remove redundant environment variables ready for next workflow
rm(list=setdiff(ls(), c('NONcancer', 'version', '%!in%', 'labkey_to_df')))

