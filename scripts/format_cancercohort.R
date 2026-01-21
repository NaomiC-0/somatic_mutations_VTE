# Script for formatting cancer cohort
# Naomi Cornish
# 19 May 2025

# Config R #######
# First run the config.R script to set up libraries and the labkey_sql function
setwd(Sys.getenv('wd'))
source('./scripts/config421.R')

# Data dictionary available here:  https://re-docs.genomicsengland.co.uk/release19/

# Query the participant_summary table
table_queried <- "participant_summary"
columns_selected <- c("participant_id", "genetically_inferred_ancestry_nothr", "programme", "date_of_consent")
# select 100K cancer programme participants
cancer_prog  <- labkey_select() %>% filter(programme == 'Cancer') %>%
	# remove duplicate records
	distinct

# confirm relevant consent status for active research ####
consent_sql <- paste("SELECT
    participant_id, programme_consent_status
    FROM participant
    WHERE programme_consent_status = 'Consenting'",
                     sep="")
active_consent <- labkey_to_df(consent_sql, version, 100000)

cancer_prog <- cancer_prog %>% filter(participant_id %in% active_consent$participant_id)

# Cancer_analysis table ####
# GEL recommend using the cancer_analysis file to identify participants where WGS has passed internal QC
table_queried <- "cancer_analysis" 
columns_selected <- c("participant_id", "tumour_sample_platekey", "germline_sample_platekey", 
                      "disease_type", "disease_sub_type", "sex", "histology_coded", "histology_coded_desc",
                      "study_name", "av_tum_date_difference", "apc_date_difference", "match_rank", "karyotype_sex",
                      "year_of_birth", "tumour_clinical_sample_time", "preparation_method", "tumour_type", "tissue_source",
                      "library_type", "at_drop", "chimeric_percentage", "average_fragment_size",
                      "tumour_sample_cross_contamination_percentage", "germline_sample_cross_contamination_percentage",
                      "tumour_purity")
cancer_analysis <- labkey_select() 

# Format this table
cancer_analysis_qc <- cancer_analysis %>%
  # tumour_sample_cross_contamination_percentage is coded as class character - group into 2 categories
  mutate(tumour_sample_cross_contamination_percentage = 
		case_when(tumour_sample_cross_contamination_percentage == '<1' ~ '<1%',
                        as.numeric(tumour_sample_cross_contamination_percentage) >= 1 & 
                        as.numeric(tumour_sample_cross_contamination_percentage) < 5 ~ '1-5%',
                        TRUE ~ NA_character_)) %>%
	# flag samples collected/stored before recruitment opened (retain only prospectively collected samples)
  # For details of recruitment process and timeline: 
  # see 'sample handling guidance v4.0' https://www.genomicsengland.co.uk/initiatives/100000-genomes-project/documentation
  mutate(prospective_sample = ifelse(tumour_clinical_sample_time < dmy('01/01/2015'), F, T)) %>%
  # **************DUPLICATES **************
  # first group by participant
  group_by(participant_id) %>%
  # for participants with multiple tumour samples deposited, select the initial/first tumour submitted to GEL.
  # if multiple tumours submitted on same initial date allow duplicates/ties
  slice_min(order_by = tumour_clinical_sample_time, with_ties =T) %>%
  # for patient who still have duplicates prioritise samples from primary tumours over metastatic site biopsy
  # convert this column to a factor first
  mutate(tumour_type = factor(tumour_type, 
                        levels=c('PRIMARY', 'RECURRENCE_OF_PRIMARY_TUMOUR', 'METASTASES'))) %>%
  # now select primary tumours
  slice_min(order_by = tumour_type, with_ties = T) %>%
  # next prioritise by tumour purity
  slice_max(order_by = tumour_purity, with_ties = T) %>%
  # next prioritise the entry where there is a matching NCRAS record
  # this will be important for determining the date of diagnosis, treatment etc
  # this value is  shown in the column match_rank: 
  # see https://re-docs.genomicsengland.co.uk/cancer_analysis_histology/#code-lists for details
  # if they have multiple entries with the same match_rank keep both for now
  slice_min(order_by = match_rank, with_ties = T) %>%
  # if duplicates remain prioritise FF (fresh frozen) samples
  mutate(preparation_method = factor(preparation_method, 
                      levels=c('FF', 'ASPIRATE', 'EDTA', 'CD128_SORTED_CELLS', 'FFPE'))) %>%
  slice_min(order_by = preparation_method, with_ties = T) %>% 
  # if duplicates still remain then select samples where the average_fragment_size is closest to the 
  # expected median value for good quality FF samples (as per the GEL data-dictionary)
  slice_min(abs(average_fragment_size - 490), with_ties = T) %>% ungroup()
  distinct

# check no duplicates following this  
cancer_analysis_qc %>% select(participant_id) %>% anyDuplicated()
rm(cancer_analysis)

# Now join the cancer_analysis_qc table to the overall GEL cancer programme
cancer <- left_join(cancer_prog, cancer_analysis_qc) %>%
# create a failed_minimum_qc flag for patients who appear in cancer_participant but not cancer_analysis
  mutate(fail_min_qc = ifelse(participant_id %!in% cancer_analysis_qc$participant_id, T, F))

# NCRAS data #####
# The National Cancer Registration and Analysis Service (NCRAS) contains information on stage and diagnosis dates.
# NB tumour_pseudo_id are assigned to participants by NCRAS and do not link to the tumour_ids assigned by GeL
# Some participants have multiple tumour samples submitted. 
# https://re-docs.genomicsengland.co.uk/cancer_analysis_histology/#cancer-analysis-histology-and-tgca-study
# this page shows how to link the NCRAS data with the cancer_analysis table 
# It is good practice to match diagnosisdatebest from av_tumour with tumour_clinicial_sample_time from cancer_analysis when merging data.
# av_tumour is the NCRAS table (prefix 'av.')

# Retrieve all the NCRAS records and then match them with correct GEL entry
table_queried <- "av_tumour" 
columns_selected <- c("participant_id", "anon_tumour_id", "age", "diagnosisdatebest", 
                      "site_coded_3char","site_coded_desc", "behaviour_coded", 
                      "behaviour_coded_desc")
av_tumour <- labkey_select() %>% 
  # according to the data dictionary, in av_tumour 'age' actually means age at diagnosis. Rename this field for clarity
  dplyr::rename('age_cancerdiag' = 'age')

# identify participants with no NCRAS record available
missing_NCRAS <- cancer %>% filter(participant_id %!in% av_tumour$participant_id)
count(missing_NCRAS)

# Merge NCRAS data with cancer_analysis data
# first retain only the participants who appear in cancer_analysis (fail_min_qc==F)
passed_gel_qc <- cancer %>% filter(cancer$fail_min_qc == F) %>% select(participant_id) %>% unlist
av_tumour_qc <- av_tumour %>% filter(participant_id %in% passed_gel_qc) %>%
  # also retain only ICD10 'C' code cancers (i.e. malignant cancers)
  filter(str_detect(site_coded_3char, 'C'))
rm(av_tumour, passed_gel_qc)

# Merging instructions as per this page:
# https://re-docs.genomicsengland.co.uk/cancer_analysis_histology/#cancer-analysis-histology-and-tgca-study
# first join av_tumour to cancer_analysis based on participant_id
cancer_ncras <- left_join(cancer, av_tumour_qc, by='participant_id') %>%
  # Create a new column named ncras_date_diff representing the absolute difference in days between diagnosisdatebest and tumour_clinical_sample_time. 
  # Replace any null values with a placeholder value of 100000 days.
  mutate(ncras_date_diff = abs(difftime(diagnosisdatebest, tumour_clinical_sample_time, units = "days")),
         ncras_date_diff = ifelse(is.na(ncras_date_diff), 100000, as.numeric(ncras_date_diff))) %>%
  # For each participant, Select the rows in the joined table with the minimum av_tum.date_difference
  group_by(participant_id) %>%
  slice_min(ncras_date_diff) %>%
  ungroup() %>%
  # create a column which flags any participants who didn't have an NCRAS record
  mutate(missing_NCRAS_record = ifelse(participant_id %in% missing_NCRAS$participant_id, T, F))

# remove the missing_NCRAS_record now as it is no longer required
rm(missing_NCRAS)

# For patients who still have duplicate records on same date, prioritise entries by histology: 3>5>2>9>6>1>0 
# 3= malignant, 5 = microinvasive, 2 = in situ, 9 = malignant, uncertain if primary or metastatic, 6= Metastatic/secondary, 1 = uncertain; 0 = benign
cancer_ncras_qc <- cancer_ncras %>% 
  # convert the behaviour column to a factor in order  of priority histology
  mutate(behaviour_coded = factor(behaviour_coded, levels = c(3, 5, 2, 9, 6, 1, 0), ordered = TRUE)) %>%
  group_by(participant_id) %>%
  # for patients with duplicate records prioritise records relating to malignant histology
 slice_min(behaviour_coded, with_ties = TRUE) %>%
  # where duplicates still persist, select unique records allowing for some minor differences in the site_coded_desc column
  # providing that the site_coded_3char is identical (this is the ICD10 code), 
  distinct(select(., -c(site_coded_desc, anon_tumour_id)), .keep_all = T) %>%
  ungroup()

# Group the cancers using the ICD10 version 19 code level 2 mappings 
# Use this to check whether the NCRAS ICD10 code is congruent with the GEL disease_type
cancer_site_fn <- function(df) {
  new_df <- df %>% 
    # first create a new column where the ICD10 C codes are converted to numeric
    mutate(numeric_code = as.numeric(gsub('C', '', site_coded_3char)),
           level2mapping = 
             case_when(numeric_code <= 14 ~ 'Head and neck',
                       numeric_code >= 30 & numeric_code <= 32 ~ 'Head and neck',
                       numeric_code >= 15 & numeric_code <= 26 ~ 'Digestive Organs',
                       numeric_code >= 33 & numeric_code <= 39 ~ 'Respiratory and Intrathoracic Organs',
                       numeric_code >= 40 & numeric_code <= 41 ~ 'Bone & Articular Cartilage and Mesothelial & Soft Tissue',
                       numeric_code >= 45 & numeric_code <= 49 ~ 'Bone & Articular Cartilage and Mesothelial & Soft Tissue',
                       numeric_code >= 43 & numeric_code <= 44 ~ 'Skin',
                       numeric_code == 50 ~ 'Breast',
                       numeric_code >= 51 & numeric_code <= 58 ~ 'Female Genital Organs',
                       numeric_code >= 60 & numeric_code <= 63 ~ 'Male Genital Organs',
                       numeric_code >= 64 & numeric_code <= 68 ~ 'Urinary Tract',
                       numeric_code >= 69 & numeric_code <= 72 ~ 'Eye, Brain and Other Parts of Central Nervous System',
                       numeric_code >= 75.1 & numeric_code <= 75.3 ~ 'Eye, Brain and Other Parts of Central Nervous System',
                       numeric_code >= 73 & numeric_code <= 75.0 ~ 'Thyroid and Other Endocrine Glands',
                       numeric_code >= 75.4 & numeric_code <= 75.9 ~ 'Thyroid and Other Endocrine Glands',
                       numeric_code >= 76 & numeric_code <= 80 ~ 'Ill-Defined, Secondary and Unspecified Sites',
                       numeric_code >= 81 & numeric_code <= 96 ~ 'Lymphoid, Haematopoietic and Related Tissue',
        # C42 doesn't appear in ICD10verision 19 so info on this code obtained from the european network of cancer registries www.encr.eu
                       numeric_code == 42 ~ 'Lymphoid, Haematopoietic and Related Tissue',  
                       numeric_code == 97 ~ 'Independent (Primary) Multiple Sites',
                       TRUE ~ NA_character_ 
                       ))
  return(new_df)}



# apply the function to the cancer_ncras_qc table
cancer_ncras_qc <- cancer_site_fn(cancer_ncras_qc)
# identify whether the GEL disease_type is congruent with the level2mapping from the ICD10 code
# the categories which do not fit the level2 mappings are CHILDHOOD and OTHER
# these will be excluded as will not be able to confirm congruence between GEL and NCRAS
cancer_ncras_qc %>% filter(disease_type == 'CHILDHOOD') %>% select(disease_sub_type) %>% table
cancer_ncras_qc %>% filter(disease_type == 'OTHER') %>% select(disease_sub_type) %>% table

# create a new column called congruent
cancer_ncras_qc <- cancer_ncras_qc %>% 
  mutate(congruent_diag = case_when(level2mapping == 'Head and neck' & disease_type == 'ORAL_OROPHARYNGEAL' ~ 'yes',
                                    level2mapping == 'Head and neck' & disease_type == 'SINONASAL' ~ 'yes',
                                    level2mapping == 'Head and neck' & disease_type == 'NASOPHARYNGEAL' ~ 'yes',
                                    level2mapping == 'Digestive Organs' & disease_type == 'UPPER_GASTROINTESTINAL' ~ 'yes',
                                    level2mapping == 'Digestive Organs' & disease_type == 'COLORECTAL' ~ 'yes',
                                    level2mapping == 'Digestive Organs' & disease_type == 'HEPATOPANCREATOBILIARY' ~ 'yes',
                                    level2mapping == 'Respiratory and Intrathoracic Organs' & disease_type == 'LUNG' ~ 'yes',
                                    level2mapping == 'Bone & Articular Cartilage and Mesothelial & Soft Tissue' & disease_type == 'SARCOMA' ~ 'yes',
                                    level2mapping == 'Skin' & disease_type == 'MALIGNANT_MELANOMA' ~ 'yes',
                                    level2mapping == 'Breast' & disease_type == 'BREAST' ~ 'yes',
                                    level2mapping == 'Female Genital Organs' & disease_type == 'OVARIAN' ~ 'yes',
                                    level2mapping == 'Female Genital Organs' & disease_type == 'ENDOMETRIAL_CARCINOMA' ~ 'yes',
                                    level2mapping == 'Male Genital Organs' & disease_type == 'PROSTATE' ~ 'yes',
                                    level2mapping == 'Male Genital Organs' & disease_type == 'TESTICULAR_GERM_CELL_TUMOURS' ~ 'yes',
                                    level2mapping == 'Urinary Tract' & disease_type == 'RENAL' ~ 'yes',
                                    level2mapping == 'Urinary Tract' & disease_type == 'BLADDER' ~ 'yes',
                                    level2mapping == 'Eye, Brain and Other Parts of Central Nervous System' & disease_type == 'ADULT_GLIOMA' ~ 'yes',
                                    level2mapping == 'Thyroid and Other Endocrine Glands' & disease_type == 'ENDOCRINE' ~ 'yes',
                                    level2mapping == 'Ill-Defined, Secondary and Unspecified Sites' & disease_type == 'CARCINOMA_OF_UNKNOWN_PRIMARY' ~'yes',
                                    level2mapping == 'Lymphoid, Haematopoietic and Related Tissue' & disease_type == 'HAEMONC' ~ 'yes',
                                    TRUE ~ NA_character_
                                    )) %>%
  # recode the congruent_diag column as a T/F factor where NA = F, otherwise TRUE.
  # note this means both NA records and congruent diagnoses will be coded as FALSE
  mutate(congruent_diag = ifelse(is.na(congruent_diag), FALSE,
                                 ifelse(congruent_diag == 'yes', TRUE, NA)))

# For the few remaining duplicate entries, prioritise rows where congruent_diag == T
cancer_ncras_qc_unique <- cancer_ncras_qc %>% 
  group_by(participant_id) %>%  # Group by participant_id
  filter(if (n() > 1) congruent_diag == TRUE else TRUE) %>%  # If participant is duplicated, keep where congruent_diag == TRUE, otherwise keep all
# flag any remaining (unresolved) duplicates for removal
  mutate(unresolved_duplicate = ifelse(n() > 1, TRUE, FALSE)) %>%  
  slice(1) %>%  # Keep only the first occurrence of any duplicate participant (they will be excluded at next step)
  ungroup()  # Ungroup to return to the original structure

# confirm no duplicates
cancer_ncras_qc_unique %>% select(participant_id) %>% anyDuplicated

# create list of GEL disease_types
gel_disease_grps <- cancer_ncras_qc_unique %>% select(disease_type) %>% distinct %>% unlist
# cycle through each disease type and where I have flagged that the diagnosis is congruent, check what the ICD10 codes are in that group.
congruent_check <- list()
for (i in gel_disease_grps) {
  congruent_check[[i]] <- cancer_ncras_qc_unique %>% 
    filter(congruent_diag == T & disease_type == i) %>% 
    count(site_coded_desc) %>% add_row(site_coded_desc = 'Total', n=sum(.$n))
    }
# to print to screen
for(name in names(congruent_check)) {print(name)
                                      print(congruent_check[[name]])
                                      # total values
                                      if (!is.null(congruent_check[[name]])) { 
                                      as.data.frame(congruent_check[[name]]) %>% filter(site_coded_desc == 'Total') %>%
                                      print}}

# to look only at the total values in each category
for(name in names(congruent_check)) {print(name)
  # total values
  if (!is.null(congruent_check[[name]])) { 
  as.data.frame(congruent_check[[name]]) %>% filter(site_coded_desc == 'Total') %>% print}}


# now check which codes did not make it through the QC step
incongruent_check <- list()
for (i in gel_disease_grps) {
  incongruent_check[[i]] <- cancer_ncras_qc_unique %>% 
    filter(congruent_diag == F & disease_type == i) %>% 
    count(site_coded_desc) %>% add_row(site_coded_desc = 'Total', n=sum(.$n))}
# print to screen
for(name in names(incongruent_check)) {print(name)
  print(incongruent_check[[name]])
  # total values
  if (!is.null(incongruent_check[[name]])) { 
  as.data.frame(incongruent_check[[name]]) %>% filter(site_coded_desc == 'Total') %>% print}}

# to look only at the total values in each category
for(name in names(incongruent_check)) {print(name)
  # total values
  if (!is.null(incongruent_check[[name]])) {
  as.data.frame(incongruent_check[[name]]) %>% filter(site_coded_desc == 'Total') %>% print}}

# clear environment
rm(cancer, cancer_ncras_qc, av_tumour_qc, cancer_analysis_qc, cancer_prog, congruent_check)

#Treatment #####
# av_treatment table: contains NCRAS  treatment data
 table_queried <- "av_treatment" 
 columns_selected <- c("participant_id", "anon_tumour_id", "eventcode", "eventdesc", "eventdate", 
                      "opcs4_code", "opcs4_name", "radiocode", "radiodesc",
                      "imagingcode", "imagingdesc", "imagingsite")
 av_treatment <- labkey_select() %>%
# Note each patient appears in this table multiple times (if they received more than one treatment)
  #convert any blank values to NA
  mutate(across(c(1:ncol(.)) & where(is.character), ~na_if(.,""))) 

# This table lists many different forms of treatment including different forms of SACT, active monitoring, radiotherapy etc
table(av_treatment$eventdesc, useNA = 'always')
# For the purpose of co-variate adjustment categorise these into either
# SACT = drug treatment including chemo / immunotherapy / hormonal therapy etc
# Surgery = any form of surgery
# NA for anything else

# identify first ever date of SACT. Note info on subsequent SACT cycles will be obtained from SACT database
first_sact <- av_treatment %>%
  # select just the participant_id, eventcode, eventdesc and eventdate columns 
  # for now the columns with more detail on the specific 'event' (treatment) are not used e.g. opcs4_code, radiocode, imagingcode etc
  select(participant_id, eventcode, eventdesc, eventdate) %>%
  # remove duplicates
  distinct %>%
  # remove the rows which relate only to imaging and surgery
  filter(!eventdesc %in% c('Imaging', 'Radiosurgery', 'surgery')) %>%
  # create a new column 'treatment_category' which groups the treatments into one of the above categories
  mutate(treatment_category = case_when(
                                        str_detect(eventdesc,'Anti-cancer drug regimen') ~ 'sact',
                                        str_detect(eventdesc,'Biological Therapies') ~ 'sact',
                                        str_detect(eventdesc,'Chemoradiotherapy') ~ 'sact',
                                        str_detect(eventdesc,'Chemotherapy') ~ 'sact',
                                        str_detect(eventdesc,'Hormone Therapy') ~ 'sact',
                                        str_detect(eventdesc,'Immunotherapy') ~ 'sact',
                                        str_detect(eventdesc,'Stem Cell Therapy') ~ 'sact',
                                        TRUE ~ NA_character_
                                        )) %>%
  # filter just for sact
  filter(treatment_category == 'sact') %>%
# For each participant, group duplicate entries by treatment category and select the first recorded entry
# Note will obtain information on subsequent SACT cycles from the SACT database lower down
# If there are duplicate entries for the same treatment category on the same date, retain the first row only
  group_by(participant_id, treatment_category) %>% slice_min(eventdate, with_ties = F) %>%  
  ungroup() %>% 
    # remove the eventcode and eventdesc columns as these are no longer required
    select(-c(eventcode, eventdesc)) %>%
    # Now reformat the dataframe so that each treatment category becomes a column with the associated eventdate as its value
  pivot_wider(names_from = treatment_category, values_from = eventdate) %>%
# rename columns for clarity
  rename(
         'date_first_sact' = 'sact',
	)

#check no duplicated participants
first_sact %>% select(participant_id) %>% anyDuplicated

# now merge treat_qc(the formatted treatment file) with the cancer_qcd file
cancer_ncras_qc_unique_first_sact <- left_join(cancer_ncras_qc_unique, first_sact, by = 'participant_id')

# select only surgical events
surgery <- av_treatment %>%
  # select just the participant_id, eventcode, eventdesc and eventdate columns 
  # for now the columns with more detail on the specific 'event' (treatment) are not used e.g. opcs4_code, radiocode, imagingcode etc
  select(participant_id, eventdesc, eventdate) %>%
  # remove duplicates
  distinct %>%
  # filter only for surgery
  filter(str_detect(eventdesc, 'Surgery'))

# To figure out when surgery occurred in relation to study entry:
# study_entry here is tumour_clinical_sample_time
# merge the surgery dates with overall cancer cohort using participant id
cancer_ncras_qc_unique_first_sact_surgery <- left_join(cancer_ncras_qc_unique_first_sact, surgery, by = join_by(participant_id)) %>%
  # note participants who received treatments on different dates will now be duplicated across rows
  # for each row (i.e. each treatment and date; calculate the time from study entry (tumour biopsy date) to sact)
  mutate(time_to_surgery = as.numeric(difftime(eventdate, tumour_clinical_sample_time, units = 'days'))) %>%
  # remove redundant columns
  select(-c(eventdate, eventdesc))

# for each participant select the surgery date which occurred closest to study entry)
cancer_ncras_qc_unique_first_sact_surgery_unique <- cancer_ncras_qc_unique_first_sact_surgery %>%
  # group by each participant  
  group_by(participant_id) %>%
  filter(
    # keep any participants where time_to_surgery is NA 
    if (all(is.na(time_to_surgery))) {
      TRUE 
    } else if (any(time_to_surgery >= -42)) {
      # if time_to_surgery is not NA
      # select the row with the minimum value of time_to_surgery, providing time_to_surgery is >= -42 (within 6 weeks study entry)
      # in other words the first surgery which occurred after -6 weeks prior to study entry
      time_to_surgery == min(time_to_surgery[time_to_surgery >= -42], na.rm = TRUE)
    } else {
      #  If time_to_surgery is < -42, and there is no row for that participant where time_to_surgery is >= -42, then select the maximum value
      # i.e. the surgery before study entry (closest to study entry)
      time_to_surgery == max(time_to_surgery, na.rm = TRUE)
    }
  ) %>%
  ungroup()  

# Stage of cancer ####
table_queried <- "cancer_staging_consolidated" 
columns_selected <- c("participant_id", "tumour_sample_platekey", "stage_best", "stage_best_system")
stage <- labkey_select()
# note stage_best_system is the system used to record best registry stage at diagnosis
# stage_best is the best 'registry' stage at diagnosis of the tumour

# merge with the main df on participant_id and tumour_sample_platekey
# this ensures that the staging entry matches the correct cancer episode for each participant
stage_cancer_ncras_treatment <- left_join(cancer_ncras_qc_unique_first_sact_surgery_unique, stage, by = c('participant_id', 'tumour_sample_platekey')) 
rm(stage)

# Metastatic disease ####   
# This table shows whether patients were recorded as having metastatic disease at diagnosis
# It is complementary to the consolidated staging table and can be used to supplement and confirm information from the staging table
table_queried <- "cancer_participant_tumour_metastatic_site" 
columns_selected <- c("participant_id", "cancer_participant_tumour_sk", "metastatic_site", "valid_metastatic_site_enumeration")
mets <- labkey_select()
# https://re-docs.genomicsengland.co.uk/release3/#cancer-participant-tumour-metastatic-site
# for my analysis only the presence or absence of metastases is relevant, not the site of metastases
# therefore I simply need to identify patients who had metastatic disease confirmed
met_confirmed <- mets %>% filter(valid_metastatic_site_enumeration == T) %>% select(participant_id) %>%
  distinct %>% unlist

# Now add a column to the main df identifying whether the participant had metastatic disease confirmed
stage_cancer_ncras_treatment <- stage_cancer_ncras_treatment %>% 
  mutate(metastatic_disease = case_when(participant_id %in% met_confirmed ~ 'yes',
                                        TRUE ~ NA_character_))

# Bin_stages #####
# Edit the staging columns so that there are fewer categories
# look at the different values for staging
table(stage_cancer_ncras_treatment$stage_best, useNA = 'always')
table(stage_cancer_ncras_treatment$stage_best_system, useNA = 'always')
# The different staging systems present in this table include:
# UICC5,6,7,8; AJCC; FIGO system,and some others. 
# Most systems can be harmonised onto a stage 1-4 scale (https://www.cancerresearchuk.org/about-cancer/what-is-cancer/stages-of-cancer#what)
# Some haem-onc staging systems: e.g. Binet and ISS and some staging systems for other rarer cancers are more difficult to harmonise. 

# https://www.cancerresearchuk.org/about-cancer/what-is-cancer/stages-of-cancer#what
# stage 0 = in situ neoplasms
# stage 1 = early (contained within organ)
# stage 2 /3 = locally invasive
# stage 4 = distant metastases
# code '?' insufficient information
# code 'U' = unstageable

# Create a list of all possible combinations for each type of staging system
staging_systems <- stage_cancer_ncras_treatment %>% distinct(stage_best_system) %>% unlist
print(staging_systems)
stages <- list()
for (i in staging_systems) {stages[[i]] <- stage_cancer_ncras_treatment %>% filter(str_detect(stage_best_system, i)) %>% select(stage_best) %>% table(., useNA = 'always')
                              }
names(stages) <- staging_systems
# to view
print(stages)

stage_cancer_ncras_treatment <- stage_cancer_ncras_treatment %>% 
  # first extract the first character from the stage_best column
  mutate(stage_temp = substring(stage_best, 1,1)) %>%
  # exclude any codes from the ISS staging system as this is the only numeric staging system which  doesn't harmonise with the numeric staging above
  mutate(stage_temp2 = ifelse(stage_best_system == 'ISS', NA, stage_temp)) %>%
  # exclude any codes which are not 0,1,2,3 or 4
  mutate(stage_numeric = case_when(stage_temp2 == '0' ~ '0',
                                   stage_temp2 == '1' ~ '1',
                                   stage_temp2 == '2' ~ '2',
                                   stage_temp2 == '3' ~ '3',
                                   stage_temp2 == '4' ~ '4',
                                   TRUE ~ NA_character_)) %>% 
  mutate(stage_numeric = as.integer(stage_numeric)) %>%
  # if stage_numeric is NA but the participant has a metastatic disease flag then populate this with stage 4
  mutate(stage_numeric = ifelse(is.na(stage_numeric) & metastatic_disease == 'yes', 4, stage_numeric)) %>%
  # remove the temporary columns
  select(-c(stage_temp, stage_temp2, stage_numeric))

# Principle components ####

table_queried <- "aggregate_gvcf_sample_stats" 
columns_selected <- c("participant_id", "pc1", "pc2", "pc3", "pc4",
                      "pred_african_ancestries", "pred_south_asian_ancestries",
                      "pred_east_asian_ancestries", "pred_european_ancestries",
                      "pred_american_ancestries")
PCs <- labkey_select()
 # select only PCs 1 - 4 as these explain most of the germline genetic variation in the cohort
# merge with the master df
pcs_stage_cancer_ncras_treatment <- stage_cancer_ncras_treatment %>% 
  left_join(., PCs, by = 'participant_id') %>%
  # ensure blank values are NA
  mutate(across(c(1:ncol(.)) & where(is.character), ~na_if(.,"")))

# Date death ####
# from ONS data
# pull death dates from labkey: participant summary table 
table_queried <- "participant_summary" 
columns_selected <- c("participant_id", "death_date")
death <- labkey_select() %>% distinct

#check for duplicates and resolve any data discrepancies if present
death %>% group_by(participant_id) %>% filter(n()>1)

pcs_stage_cancer_ncras_treatment <- death %>% 
  left_join(pcs_stage_cancer_ncras_treatment, ., by = 'participant_id')

#Save #####
fwrite(pcs_stage_cancer_ncras_treatment, './data/pcs_imd_stage_cancer_ncras_treatment.csv', row.names = F)


# END OF SCRIPT ####





