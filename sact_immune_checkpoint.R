# Format SACT data
# N Cornish

#Here I specifically identify patients who received immune checkpoint inhibition
# Config ####
setwd(Sys.getenv('wd'))
source('./scripts/config421.R')
# Load formatted analysis cohort ####
ac <- readRDS('./data/analysis_cohort_formatted.rds')

# Query SACT database ####
table_queried <- "sact"
columns_selected <- c("participant_id", "anon_tumour_id", "primary_diagnosis", 
                      "start_date_of_regimen", "drug_group")

sact <- labkey_select()

fwrite(sact, './data/sact.csv', row.names = F)

# to read back in
sact <- fread('./data/sact.csv') %>%  # ensure blank character values are coded NA    
  mutate(across(c(1:ncol(.)) & where(is.character), ~na_if(.,"")))

# identify immune check point inhibitors
#obtained from here: https://www.cancerresearchuk.org/about-cancer/treatment/targeted-cancer-drugs-immunotherapy/checkpoint-inhibitors
# PD-1: nivolumab (Opdivo); pembrolizumab (Keytruda); cemiplimab (Libtayo)
# CTLA-4 inhibitors: Ipilimumab (Yervoy)
# PD-L1 inhibitors: atezolizumab; avelumab; durvalumab
# LAG-3: Relatlimab 

# CDK inhibitors: abemaciclib, palbociclib, ribociclib


immune_checkpoint <- sact %>% 
  filter(str_detect(drug_group, 'NIVOLUMAB')|str_detect(drug_group, 'OPDIVO')|
         str_detect(drug_group, 'PEMBROLIZUMAB')|str_detect(drug_group, 'KEYTRUDA')|
           str_detect(drug_group, 'CEMIPLIMAB')|str_detect(drug_group, 'LIBTAYO')|
           str_detect(drug_group, 'IPILIMUMAB')|str_detect(drug_group, 'YERVOY')|
           str_detect(drug_group, 'ATEZOLIZUMAB')|str_detect(drug_group, 'AVELUMAB')|
           str_detect(drug_group, 'DURVALUMAB')|str_detect(drug_group, 'RELATLIMAB')
         ) %>% 
  distinct(drug_group) %>% pull %>% as.data.frame %>% mutate(category = 'immune_checkpoint') %>% rename('drug' = '.')

cdk_inhibitors <- sact %>% filter(str_detect(drug_group, 'ABEMACICLI')|
                                    str_detect(drug_group, 'PALBOCICLIB')|
                                   str_detect(drug_group, 'RIBOCICLI'))%>% 
  distinct(drug_group) %>% pull %>% as.data.frame %>% mutate(category = 'CDK_inhibitor') %>% rename('drug' = '.')

sact_categories2 <- rbind(cdk_inhibitors, immune_checkpoint)
  
# Look at the most common drug groups
sact_grps <- sact %>% 
  # annotate the sact dataframe with the umbrella categories which I have created
  left_join(., sact_categories2, by = join_by(drug_group == drug)) %>%
  # exclude any rows that are not SACT (i.e. not in the sact_categories2 table)
  filter(drug_group %in% sact_categories2$drug) %>%
  # filter only for participants in the analysis cohort
  filter(participant_id %in% ac$participant_id) %>%
  # group by drug category and participant
  group_by(participant_id, category) %>% 
  # ensure that each participant is only represented once for each drug category, then table
  slice(1) %>% ungroup %>% select(category) %>% table %>% as.data.frame %>%
  setnames(., colnames(.), c('drug category', 'freq')) %>% arrange(desc(freq))

write.table(sact_grps, './results/sact_freqs_immunos_CDKi.txt', row.names=F, quote = F, sep = '\t')

# join the SACT df with the sact_categories2 so that only CKDi and immune checkpoint inhibitors are retained
sact_subset <- left_join(sact_categories2, sact, by = c('drug' = 'drug_group')) %>%
  relocate(c(drug, category), .after = last_col())

# time to therapy ####
# merge with ac and calculate time to CDKi and time to immune checkpoint inhibition

ac_cdk <- sact_subset %>% filter(category == 'CDK_inhibitor') %>%
  # select only required columns (participant_id, sact and drug group for clarity)
  select(participant_id, start_date_of_regimen, category) %>%
  # for participants with multiple prescriptions select earliest date
  group_by(participant_id) %>% slice_min(start_date_of_regimen) %>%
  # remove any duplicate rows
  distinct %>%
  # merge with the overall cancer cohort using participant id
  left_join(ac, ., by = join_by(participant_id)) %>%
  # note participants who received treatments on different dates will now be duplicated across rows
    # for each row (i.e. each treatment and date; calculate the time from study entry to sact)
  mutate(time_to_CDKi = as.numeric(difftime(start_date_of_regimen, study_entry, units = 'days')))

# to check
ac_cdk %>% filter(!is.na(time_to_CDKi)) %>% select(category) %>% table(., useNA = 'always')


ac_ICI <- sact_subset %>% filter(category == 'immune_checkpoint') %>%
  # select only required columns (participant_id, sact and drug group for clarity)
  select(participant_id, start_date_of_regimen, category) %>%
  # for participants with multiple prescriptions select earliest date
  group_by(participant_id) %>% slice_min(start_date_of_regimen) %>%
  # remove any duplicate rows
  distinct %>%
  # merge with the overall cancer cohort using participant id
  left_join(ac, ., by = join_by(participant_id)) %>%
  # note participants who received treatments on different dates will now be duplicated across rows
  # for each row (i.e. each treatment and date; calculate the time from study entry to sact)
  mutate(time_to_ICI = as.numeric(difftime(start_date_of_regimen, study_entry, units = 'days')))

# to check
ac_ICI %>% filter(!is.na(time_to_ICI)) %>% select(category) %>% table(., useNA = 'always')

# merge them, just retaining participant ID and time_to_ICI or time_to_CDKi

ac_ICI_CDKi <- ac_cdk %>% select(participant_id, time_to_CDKi) %>% 
  left_join(., ac_ICI %>% select(participant_id, time_to_ICI), 
            by = 'participant_id')
                                                                             

rm(ac_ICI, ac_cdk)

# save
fwrite(ac_ICI_CDKi, './data/ac_ICI_CDKi.csv', row.names = F)

