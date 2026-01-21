# Format analysis cohort
# N Cornish

# Notes ######
# Script migrated from main workflow in order to streamline it for github
# It was originally part of a longer Rmd document and I have not re-run in this new format
# Consequently there may be some bugs

source('config421.R')

# Merge into single df #######
# Load the relevant tables created in steps 1-3 and merge into single df
# output from ./format_cancer_cohort.R
cancer <- fread('./data/pcs_imd_stage_cancer_ncras_treatment.csv') %>%
  # ensure blank character values are coded NA    
  mutate(across(c(1:ncol(.)) & where(is.character), ~na_if(.,""))) %>%
  # rename the date_first_sact column to make clear that it comes from the av_treatment table rather than the sact table as there are minor discrepancies between the dates recorded in these two systems
  rename('date_first_sact_from_av_treatment' = date_first_sact)

# output from ./sact_data.R
sact_data_formatted <- fread('./data/sact_formatted_cancercohort.csv')

# merge
cancer_sact <-left_join(cancer, 
                        sact_data_formatted)

# output from ./VTE_diagnoses_cancer.R
vte_cancercohort <- fread('./data/hes_combined_VTE_cancercohort.csv')

# merge
dat  <- left_join(cancer_sact, vte_cancercohort, by = 'participant_id')

# Set status columns ####
# Set censor date
last_follow_up <- dmy('31/07/2022')

dat <- dat %>%
  mutate(study_entry = tumour_clinical_sample_time, 
         # calculate time to VTE (if applicable) from the observation start time
         time_to_vte = as.numeric(difftime(as.Date(date_vte), study_entry, units = 'days')),
         # calculate time to death (if applicable) from observation start time
         time_to_death = as.numeric(difftime(as.Date(death_date), study_entry, units = 'days')),
         # calculate time to last follow up 
         time_to_last_follow_up = as.numeric(difftime(as.Date(last_follow_up), study_entry, units = 'days')),
         # calculate administrative censoring: here set at 5 years from study entry
         time_to_admin_censor = 5*365) %>%
  # select the minimum value from the different time columns (ignoring any NA values)
  rowwise() %>%
  mutate(time_to_status = min(c(time_to_vte, time_to_death, time_to_last_follow_up,
                                time_to_admin_censor), na.rm = T)) %>% ungroup %>%
  # assign a status
  # if the vte was the first event and occurred on or after the observation start, code as vte
  mutate(
    event_vte = ifelse(time_to_status == time_to_vte & time_to_status >= 0, 1, 0),
    # if the first event was death, code as death (participants who died on same day as VTE diagnosis should be coded as VTE)
    event_death = ifelse(time_to_status == time_to_death, 1, 0),
    # if the first event was last follow up or admin censoring, code as censor
    event_censor = ifelse(time_to_status == time_to_last_follow_up|time_to_status == time_to_admin_censor, 1, 0),
    # if a vte occurred before the observation start time code as exclude
    event_exclude = ifelse(time_to_status < 0, 1, 0 )) %>%
  # Now assign an overall status
  # code 1 = vte, code 2 = death (competing risk), code 0 = censored/no VTE, NA = prior VTE (excluded)
  mutate(status = as.factor(case_when(event_censor == 1 ~0,
                                      event_vte == 1 ~ 1,
                                      event_death == 1 ~2,
                                      event_exclude == 1 ~ 9))) %>%
  mutate(status_description = case_when(status == 0 ~ 'no vte',
                                        status == 1 ~ 'vte',
                                        status == 2 ~ 'died',
                                        status == 9 ~ 'exclude, prior vte',)) %>%
  mutate(status_description = factor(status_description, levels = c('no vte', 'vte', 'died', 'exclude, prior vte'))) %>% 
  # convert the time to status to months for easier display of the data
  mutate(time_to_status_months= round(time_to_status/30.417, digits = 0))

# Exclusions #####
# Total number of particicpants in cancer programme with ongoing consent for research. 
cancer_prog_consent <-count(dat) %>% unlist

dat_qc <- dat %>% 
  # flag any samples where the tumour_cross_contamination is not <5% as failing minimum QC; otherwise keep the original value for fail_min_qc
  mutate(fail_min_qc = ifelse
         # due to the way this field was coded in GEL (mixture of character '<1' and numeric values) it has been binned and converted to a factor here
         (tumour_sample_cross_contamination_percentage != '<=1%' & tumour_sample_cross_contamination_percentage != '1-5%', T, fail_min_qc)) %>%
  # flag samples with discordant phenotypic and karyotypic sex
  mutate(concordant_sex = ifelse(sex == 'FEMALE' & karyotype_sex == 'XX', T,
                                 ifelse(sex == 'MALE' & karyotype_sex == 'XY', T, F))) %>%
  mutate(fail_min_qc = case_when(concordant_sex != T ~ T, 
                                 .default = fail_min_qc))

# Count the total number now failing QC 
fail_qc <- dat_qc %>% filter(fail_min_qc == T) %>% count %>% unlist
# number remaining after excluding failed QC samples
include1 <- dat_qc %>% filter(fail_min_qc == F)
rm(dat_qc)

# exclude participants where stored 'retrospective' tumour biopsies were used for WGS (this is defined as samples collected prior to recruitment opening in Jan 2015)
retrospec_sample <- include1 %>% filter(prospective_sample == F) %>% count %>% unlist

include2 <- include1 %>% 
  filter(prospective_sample == T)
rm(include1)

# no NCRAS linkage
no_ncras <- include2 %>% filter(missing_NCRAS_record == T) %>% count %>% unlist

# NCRAS linkage but incongruent diagnosis between GEL and NCRAS
incong <- include2 %>% filter(missing_NCRAS_record == F & congruent_diag == F) %>% count %>% unlist

# NCRAS linkage and congruent diagnosis but non malignant histology
non_malig <- include2 %>% filter(missing_NCRAS_record == F & congruent_diag == T) %>% 
  # filter for non-C code or BENIGN histology
  filter(!str_detect(site_coded_3char, 'C')|behaviour_coded_desc == 'BENIGN') %>% count %>% unlist

# inclusion set
include3 <- include2 %>% filter(missing_NCRAS_record == F & congruent_diag == T & str_detect(site_coded_3char, 'C') & behaviour_coded_desc != 'BENIGN')
rm(include2)

# participants with missing/discrepant data for critical covariates including genetic princial components
missing_covars <- include3 %>% 
  # filter for missing ancestry, diagnosis date, age, sex or for the small number who had duplicated records with data discrepancies which could not be resolved
  filter(is.na(genetically_inferred_ancestry_nothr)|is.na(diagnosisdatebest)|is.na(sex)|is.na(age_cancerdiag)|unresolved_duplicate == T|death_date_conflict == T) %>%
  count %>% unlist 


include4 <- include3 %>%  filter(!is.na(genetically_inferred_ancestry_nothr) & !is.na(diagnosisdatebest) & !is.na(sex) & !is.na(age_cancerdiag)) %>%
  filter(unresolved_duplicate == F)
rm(include3)

# VTE prior to study entry
prior_vte <- include4 %>% filter(status == 9) %>% count %>% unlist

# of the people who were excluded due to a prior history of VTE, identify when these events occurred

# people with a VTE before cancer diagnosis 
vte_prior_diagnosis <- include4 %>% filter(date_vte <= diagnosisdatebest) %>% filter(date_vte < study_entry) %>% count %>% unlist

# people with a VTE between cancer diagnosis and study entry
vte_between_diag_and_study_entry <- include4 %>% filter(date_vte < study_entry) %>% filter(date_vte > diagnosisdatebest) %>% count %>% unlist

# analysis cohort
analysis_cohort1 <- include4 %>% filter(status != 9)
ac_count <- count(analysis_cohort1) %>% unlist

# Sensitivity analyses #######

#prior_anticoag_indication column flags participants with a medical indication for anticoagulation before study entry

# input file: list of anticoagulation indications: see supplementary table in manuscript
source(anticoag_indications.R)
# output file: ./data/hes_combined_anticoag_cancercohort.csv
anticoag_participants <- fread('./data/hes_combined_anticoag_cancercohort.csv')

anticoag_before_studyentry <- analysis_cohort1 %>%
  #first select relevant columns from the analysis cohort
  select(participant_id, study_entry) %>%
  # merge with the anticoagulation data
  left_join(anticoag_participants, .) %>%
  # filter for participants where anticoag date is before study entry
  filter(date_anticoag <= study_entry)

# Flag these individuals in the analysis_cohort to exclude them in a sensitivity analysis
analysis_cohort2 <- analysis_cohort1 %>% 
  mutate(prior_anticoag_indication = case_when(participant_id %in% anticoag_before_studyentry$participant_id ~ 'yes',
                                               .default = 'no'))

#new_untreated_cancer column: participants who had tumour biopsied for WGS within 42 days of initial cancer diagnosis AND did not have any SACT before study entry.
# for this cohort, time zero is set at cancer diagnosis
# participants who experienced a vte in between diagnosis and study entry are included (use include4 to create this)

ac_sensitivity_diag_equals_time0 <- include4 %>% 
  # allow anyone from include 4 who had a VTE on/after cancer diagnosis, or on/after study entry (in case tumour biopsy occurred before diagnosis) and also anyone who did not have a VTE. 
  filter(date_vte >= diagnosisdatebest|date_vte >= study_entry|is.na(date_vte)) %>%
  # people with VTE in between cancer diagnosis and study entry are coded 9 here 
  # need to edit the status column for these people: they should now be 1
    mutate(status = case_when(status == '9' ~ '1', .default = status),
         # Step 1: Calculate time from diagnosis to study entry (in same units as time_to_status)
         # This is the "delay" between diagnosis and recruitment
         diag_to_entry = as.numeric(difftime(study_entry, diagnosisdatebest, units = "days")),
         # Step 2: Calculate new time from diagnosis to status
         # Add the delay (diag_to_entry) to the original time_to_status
         # But only if diagnosis occurred before study entry (diag_to_entry > 0)
         time_from_diagnosis_to_status = if_else(
           diagnosisdatebest < study_entry,
           diag_to_entry + time_to_status,  # diagnosis before entry: add the delay
           time_to_status                    # diagnosis after/at entry: keep original
         )
  ) %>%
  # identify then new_untreated_cancer group below
  # list of participants which excludes chemo prior to enrolment; diagnosis date > 42 days (6weeks) before enrolment 
# note for this df I already created the diag_to_entry column above
dplyr::rename(time_diag_to_biopsy = diag_to_entry) %>%
  mutate(time_to_first_sact_av_treatment = (as.numeric(difftime(date_first_sact_from_av_treatment, study_entry, units = "days")))) %>%
  # if a participant had time_diag_to_biopsy > 42 days or received prior sact marks as 'no' in 'new_untreated_cancer' column
  mutate(new_untreated_cancer = case_when(time_diag_to_biopsy > 42|
                                            sact_over6weeks_before_studyentry == T|
                                            time_to_first_sact_av_treatment < 0
                                          ~ 'no', .default = 'yes'))


new_untreated_count <- ac_sensitivity_diag_equals_time0 %>% filter(new_untreated_cancer == 'yes') %>% count %>% unlist

# STROBE diagram values ######

# This code derives the numbers of participants excluded from the analysis due to the reasons outlined in the methods. See Figure 1 of manuscript.
# consort values obtained from the exclusions section above
exclusions <- data.frame('all_consenting' = cancer_prog_consent, 'fail_qc' = fail_qc, 'tumour_sample_prior_2015' = retrospec_sample, 'missing_NCRAS' = no_ncras, 'incongruent_diagnosis' = incong, 'benign' = non_malig,
                         'missing_or_discrepant_covariates' = missing_covars, 'prior_vte'= prior_vte,
                         'prior_vte_before_diagnosis' = vte_prior_diagnosis, 'prior_vte_between_diag_and_study' = vte_between_diag_and_study_entry,
                         'primary_cohort'= ac_count) %>% 
  # pivot the table so that each value appears on a single row
  pivot_longer(., cols = everything(), names_to = 'group', values_to = 'count') %>%
  # calculate percentages of the original GEL cohort
  mutate(percent = paste0('(',round((count/ cancer_prog_consent)*100, digits=1),'%)')) %>%
  mutate(`n (%)` = paste(count, percent))

write.table(exclusions, './results/STROBE_diagram_vals.tsv', row.names= F, quote = F, sep = '\t')

rm(exclusions, sensitivity_cohorts, consort_values, ac_count, non_malig, missing_covars, no_ncras, prior_vte, incong, fail_qc, retrospec_sample, staging_data_complete, new_untreated_count, prior_anticoag_indication_count, anticoag_before_studyentry, anticoag_participants)

# Format column classes ####
# format column classes (factor, date etc) 
# derives new columns where required so that the table can be fed into various different statistical models with minimal additional formatting.

ac1 <- analysis_cohort2 %>% 
  # calculate age at study entry
  mutate(age_study_entry = as.integer(format(study_entry, "%Y")) - year_of_birth)

#Group 'HEAD_NECK' tumours together and group other rare tumour types as 'OTHER'

ac2 <- ac1 %>% mutate(disease_type = 
                        case_when(disease_type == 'ORAL_OROPHARYNGEAL' ~ 'HEAD_NECK',
                                  disease_type == 'NASOPHARYNGEAL' ~ 'HEAD_NECK',
                                  disease_type == 'SINONASAL' ~ 'HEAD_NECK',
                                  disease_type == 'ENDOCRINE' ~ 'OTHER',
                                  disease_type == 'TESTICULAR_GERM_CELL_TUMOURS' ~ 'OTHER',
                                  disease_type == 'CARCINOMA_OF_UNKNOWN_PRIMARY' ~ 'OTHER',
                                  .default = disease_type))

#create a new 'binary_status' column where death is treated as a censoring event (for cox models)
ac3 <- ac2 %>% 
  mutate(binary_status = as.integer(case_when(status_description == 'no vte' ~ 0,
                                              status_description == 'died' ~ 0, 
                                              status_description == 'vte' ~ 1,)))

# Recode participants who received SACT/ surgery within 6 weeks prior to study entry as exposed to these treatments at baseline
ac4 <- ac3 %>% 
  mutate(
    #  if time_to_current_sact is the same or greater than the max period of observation for that participant convert it to NA because that is how it will be handled in the anlaysis i.e. SACT did not contribute to the exit event
    time_to_current_sact = case_when(time_to_current_sact >= time_to_status ~ NA,
                                     .default = time_to_current_sact)) %>%
  # if sact occurred prior to but within 6 weeks of study entry then change time_to_sact to zero: participant was exposed to SACT at baseline (cox model cannot handle negative time values)
  mutate(
    # note due to the way the sact data was formatted (see above), participants who had SACT > 6 weeks before study entry will be NA for time_to_current_sact
    time_to_current_sact = case_when(time_to_current_sact < 0 ~ 0,
                                     # time_to_current_sact within the observation period remain unchanged
                                     .default = time_to_current_sact)) %>%
  # format time_to_surgery column
  mutate( # if surgery occurred before -42 days or after time_to_status set as NA as surgrey did not occur in study period
    time_to_surgery = case_when(time_to_surgery < -42|
                                  time_to_surgery >= time_to_status ~ NA,
                                # otherwise use the true value
                                .default = time_to_surgery)) %>%
  # create a column identifying whether participant has had surgery during the study period
  mutate(surgery_during_study_period = case_when(!is.na(time_to_surgery) ~ T,
                                                 .default = F))

# group cancer stages into early vs advanced or unknown if staging not available
ac5 <- ac4 %>% 
  mutate(stage_grouped = case_when(stage_numeric_imputed <= 2 ~ 'early',
                                   stage_numeric_imputed >= 3 ~ 'late',
                                   is.na(stage_numeric_imputed) ~ 'unknown'))

# Identify people who had SACT during the study period
ac6 <- ac5 %>% mutate(sact_during_study = case_when(time_to_current_sact < time_to_status ~ T, .default = F))

# Select only the columns required for the analysis and ensure that columns are coded as the correct class
ac <- ac6 %>% 
  # select only relevant columns that are needed for the analysis
  select(participant_id, 
         tumour_sample_platekey,
         germline_sample_platekey,
         sex, 
         age_study_entry,
         genetically_inferred_ancestry_nothr, 
         disease_type, 
         disease_sub_type, 
         site_coded_3char,
         stage_grouped,
         pc1, pc2, pc3, pc4, 
         diagnosisdatebest, 
         study_entry,
         sact_over6weeks_before_studyentry,
         time_to_current_sact,
         time_to_surgery,
         time_diag_to_biopsy,
         sact_during_study,
         surgery_during_study_period, 
         status, status_description, binary_status,
         time_to_status, time_to_status_months, new_untreated_cancer, prior_anticoag_indication) %>%
  # drop factors that do not cocur in the data (so for example status = 9 code since this has now been excluded)
  mutate(status = droplevels(as.factor(status)),
         status_description = factor(status_description, levels = c('no vte', 'vte', 'died')),
         # convert sex, ancestry and disease_type to factors,
         # use reference factors of Female, European and 'PROSTATE' respectively
         # set reference factor of 1 for cancer stage
         sex = factor(sex, levels = c('FEMALE', 'MALE')),
         genetically_inferred_ancestry_nothr = relevel(factor (genetically_inferred_ancestry_nothr), ref= 'European'),
         disease_type = relevel(factor(disease_type), ref = 'PROSTATE'),
         stage_numeric_imputed = relevel(factor(stage_numeric_imputed), ref = '1'),
         stage_grouped= factor(stage_grouped, levels = c('early', 'late', 'unknown')))
rm(ac4, ac5, ac6)

# save as R file to retain column classes 
#saveRDS(ac, './data/analysis_cohort_formatted.rds')
