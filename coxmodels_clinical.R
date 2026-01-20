# Univariate Cox regression models looking at association between selected covariates and VTE
# N Cornish

# script for running univariate Cox regression models 

# Config ####
# load config script
setwd(Sys.getenv('wd'))
source('./scripts/config421.R')
library(survival)
library(broom)


# Load formatted analysis cohort ####
ac <- readRDS('./data/analysis_cohort_formatted.rds') %>%
  # select relevant columns
  select(participant_id, tumour_sample_platekey, # columns required to merge on genotype
         age_study_entry, age_study_entry_scaled, age_cancerdiag,
         sex, pc1, pc2, pc3, pc4, # covariates for analysis
         genetically_inferred_ancestry_nothr,
         disease_type, stage_numeric_imputed, stage_grouped, 
         study_entry, sact_over6weeks_before_studyentry, time_to_current_sact,time_to_surgery,
         surgery_during_study_period,
         binary_status, time_to_status # status and time for cox model
  ) %>%
  # recode the NA stage values as a categorical covariate
  mutate(stage_numeric_imputed = as.factor(case_when(is.na(stage_numeric_imputed) ~ 'unknown stage', 
                                           .default = stage_numeric_imputed)))

# create sdata where sact is coded as a time varying covariate
##tmerge ####
# use tmerge because sact is a time dependent covar
# first create the basic dataset for tmerge
tdata <- ac %>% 
  # add +0.5 to time_to_status only where time_to_status =0
  # this step is required becuse tmerge does not accept stop times of 0
  # see vignette by Therneau et al
  mutate(time_to_status = case_when(time_to_status == 0 ~ 0.5,
                                    .default = time_to_status ))
# now perform tmerge to split f/u times for participants who received sact after study entry
sdata <- tmerge(data1 = tdata, data2 = tdata,
                vte = event(time_to_status, binary_status), 
                current_sact = tdc(time_to_current_sact), 
                id = participant_id,
                options= list(idname="participant_id")) %>%
  # select only relevant columns
  select(participant_id, tumour_sample_platekey, age_study_entry, age_study_entry_scaled, age_cancerdiag, 
         sex, genetically_inferred_ancestry_nothr,
         disease_type, stage_numeric_imputed, stage_grouped,
         sact_over6weeks_before_studyentry, current_sact, surgery_during_study_period,
         tstart, tstop, vte)

# Univariable cox regressions ####

# create vector of covariates to query in univariate regression
covars <- sdata %>% select(age_study_entry, age_study_entry_scaled, age_cancerdiag, 
                        sex, genetically_inferred_ancestry_nothr,
                        stage_numeric_imputed, stage_grouped,
                        sact_over6weeks_before_studyentry, current_sact, surgery_during_study_period) %>% 
  names %>% unlist
# cancer type regressions run seperately later

results1 <- list()
# univariable cox models
for (i in covars) {
  formula <- as.formula(paste('Surv(tstart, tstop, vte) ~', i))
  results1[[i]] <- coxph(formula, data = sdata) %>% 
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # round pval to 2dp (while retaining scientific notation)
    mutate(p.value = formatC(p.value, format = 'e', digits=2)) %>%
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    rename('covariate' = term, 'HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(covariate, HR, l95ci, u95ci, pval)
  }

results1_df <- do.call(rbind, results1)

# Univar tumour type vs all others ######
# these cox models compare HR for VTE in each tumour type to all other tumours combined
tumour_types <- ac %>% pull(disease_type) %>% unique()
# empty list to store results
results2 <- list()
# for each GEL disease category
for(i in tumour_types) {
  # create a temp dataframe where that disease is coded as 1 and all other cancers coded as 0
  temp_data <- sdata %>% 
    mutate(cancer_type_binary = case_when(disease_type == i ~ 1, .default = 0))
  # now run the regression. Use cancer_type_binary instead of disease_type
  results2[[i]] <- coxph(Surv(tstart, tstop, vte) ~ cancer_type_binary, data = temp_data) %>% 
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # rename this value to make it clear which cancer type the HR refers to
    mutate(term = paste(i, '(ref all other cancers)')) %>%
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    rename('covariate' = term, 'HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(covariate, HR, l95ci, u95ci, pval)
}

results2_df <- do.call(rbind, results2)

# Combine all results ####

results_df <- rbind(results1_df, results2_df)%>%
  # create a new column called category for clarity
  mutate(category = case_when(str_detect(covariate, '^age') ~ 'Age',
                              str_detect(covariate, 'sex') ~ 'Sex (ref = FEMALE)',
                              str_detect(covariate, 'genetically_inferred_ancestry_nothr') ~ 'Genetic ancestry (ref = European)',
                              str_detect(covariate, '(ref all other cancers)') ~ 'Cancer type (ref all other cancers)',
                              str_detect(covariate, 'stage_numeric') ~ 'Stage of cancer (0-4, ref=1)',
                              str_detect(covariate, 'stage_grouped') ~ 'Stage of cancer (early/late/unknown, ref=early)',
                              str_detect(covariate, 'sact') ~ 'sact',
                              str_detect(covariate, 'surgery') ~ 'surgery'
                              ))%>%
  # clean the text in the covariate column to remove superfluous text since this information is now shown in category
  mutate(covariate = gsub('genetically_inferred_ancestry_nothr|sex|disease_type|_imputed|TRUE|stage_numeric', '', covariate)) %>%
  # have to remove unwanted brackets in the (ref all other cancers) use different regex
  mutate(covariate = gsub("\\s*\\([^\\)]+\\)", "", covariate)) %>%
  # convert category column to a factor so that categories do not default to alphabetical order in the plot.
  mutate(category = factor(category, levels = c('Age', 'Sex (ref = FEMALE)', 'Genetic ancestry (ref = European)',
                                                'Cancer type (ref all other cancers)', 'Stage of cancer (0-4, ref=1)',
                                                'Stage of cancer (early/late/unknown, ref=early)',
                                                'sact', 'surgery'))) %>%
  # then arrange rows first by category and then in order of HR
  arrange(category, desc(HR)) %>%
  # add a column which states the model run
  mutate(model_covars = 'univariable')

# save
write.table(results_df, './results/cox_clinical_univar.tsv', quote = F, row.names = F, col.names = T, sep = "\t")
rm(results1, results2, results_df, results1_df, results2_df)



