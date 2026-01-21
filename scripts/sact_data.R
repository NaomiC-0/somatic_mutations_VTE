# Format SACT data
# N Cornish
# 25th June 2025

# Config ####
setwd(Sys.getenv('wd'))
source('./scripts/config421.R')

# Load cancer cohort (created in format_cancer_cohort.R script) #####
# study_entry here is tumour_clinical_sample_time
cancer <- fread('./data/pcs_imd_stage_cancer_ncras_treatment.csv') %>%
  # ensure blank character values are coded NA    
  mutate(across(c(1:ncol(.)) & where(is.character), ~na_if(.,""))) %>%
  # select just key columns
  select(participant_id, disease_type, disease_sub_type, tumour_clinical_sample_time, anon_tumour_id, 
                         diagnosisdatebest, date_first_sact
                         ) %>%
  # rename columns for clarity
  rename('study_entry' = tumour_clinical_sample_time, 
         'date_first_sact_from_av_treatment' = date_first_sact) # to distinguish it from the dates in the SACT database
# note the SACT database only started in 2012 whereas the av_treatment table contains records of chemo received prior to this date

# Query SACT database ####
table_queried <- "sact"
columns_selected <- c("participant_id", "anon_tumour_id", "primary_diagnosis", 
                      "start_date_of_regimen", "drug_group")

sact <- labkey_select()

fwrite(sact, './data/sact.csv', row.names = F)

# to read back in
sact <- fread('./data/sact.csv') %>%  # ensure blank character values are coded NA    
  mutate(across(c(1:ncol(.)) & where(is.character), ~na_if(.,"")))

# SACT dates #####
# Note in the format_cancer_cohort.R script, I have recorded the first_ever date of chemo from the av_treatment table
# this is shown in the column date_first_sact_from_av_treatment
# for the sake of including chemo as a time dependent covariate, identify whether a participant either:
# a) Never had chemo,
# b) Received chemo more than 6 weeks  prior to study enrolment 
# b) Had chemo within 6 weeks of study entry: for this analysis classified as being on chemo at baseline
# c) Had chemo following study entry: in which case chemo is included as a time varying covariate

# first merge the sact data with the cancer cohort dataframe
cohort_sact <- sact %>% 
  # exclude any values which are 'NOT CHEMO', 'NOT MATCHED' or 'NA'
  filter(drug_group != 'NOT CHEMO' & drug_group != 'NOT MATCHED' & !is.na(drug_group)) %>%
  # select only required columns (participant_id, sact and drug group for clarity)
  select(participant_id, start_date_of_regimen, drug_group) %>%
  # remove any duplicate rows
  distinct %>%
  # merge with the overall cancer cohort using participant id
  left_join(cancer, ., by = join_by(participant_id)) %>%
  # note participants who received treatments on different dates will now be duplicated across rows
    # for each row (i.e. each treatment and date; calculate the time from study entry to sact)
  mutate(time_to_sact = as.numeric(difftime(start_date_of_regimen, study_entry, units = 'days')))

# identify anyone with a history of SACT more than 6 weeks (42 days) prior to study entry 
# (based on EITHER the start_date_of_regimen columns or the date_first_sact_from_av_treatment)
sact_over6weeks_before_studyentry <- cohort_sact %>% 
  # ensure study entry is handled as a date first
  mutate(study_entry = as.Date(study_entry)) %>%
  # filter any participant where the date first sact is > 42 days before study entry (either in SACT or NCRAS)
  filter(date_first_sact_from_av_treatment < (study_entry - 42) |start_date_of_regimen < (study_entry - 42)) %>%
  pull(participant_id) %>% unique

# Now add a column identifying whether people had sact before study entry
cohort_sact_unique <- cohort_sact %>%
  mutate(sact_over6weeks_before_studyentry = ifelse(participant_id %in% sact_over6weeks_before_studyentry, T, F)) %>%
# next create a 'time_to_current_sact' column which is NA if SACT only occurred >6 weeks before study entry
  mutate(time_to_current_sact = ifelse(time_to_sact >= -42, time_to_sact, NA)) %>%
  # for participants who have multiple SACT dates occuring AFTER study entry, select earliest date
  # group by each participant  
  group_by(participant_id) %>%
  filter(
    # keep any participants where time_to_current_sact is NA  
    if (all(is.na(time_to_current_sact))) {
      TRUE 
    } else {
    # if time_to_current_sact is not NA
    # select the row with the minimum value (i.e. first sact occurring after day -42 in relation to study entry)
      time_to_current_sact == min(time_to_current_sact, na.rm = TRUE)
    })%>%
  ungroup()  %>%
   # select relevant columns and remove duplicate rows 
   # note some participants have multiple different drugs prescribed on same date, but the specific regimen is not used as a covar
    select(participant_id, study_entry, sact_over6weeks_before_studyentry, time_to_current_sact) %>%
  distinct

# save
fwrite(cohort_sact_unique, './data/sact_formatted_cancercohort.csv', row.names = F)

