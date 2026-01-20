# Cox regression models for signature analyses
# N Cornish

# Config ####
# load config script
setwd(Sys.getenv("wd"))
source('./scripts/config421.R')
library(survival)
library(broom)

# Load formatted analysis cohort ####
ac <- readRDS('./data/analysis_cohort_formatted.rds')

# load the file with the TMBs and mutational signatures 
tmb_ac_continuous <- fread('./data/tmb_and_mutsig_ac.csv') %>% 
  # select only the columns relating to signatures
  select(tumour_sample_platekey, starts_with('signature'))
# convert signatures (continuous numerical with high proprtion NA/zero vals to a binary val)
  # if NA then code as zero, otherwise code as 1
tmb_ac <- tmb_ac_continuous %>% mutate(across(starts_with('signature'), ~ifelse(is.na(.), 0,1)))

# create sdata where sact is coded as a time varying covariate
##tmerge 
# use tmerge because sact is a time dependent covar
# first create the basic dataset for tmerge
tdata <- ac %>% 
  # merge with the signatures data
  left_join(., tmb_ac, by = 'tumour_sample_platekey') %>%
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
         sex, pc1, pc2, pc3, pc4,
         disease_type, stage_numeric_imputed, stage_grouped,
         sact_over6weeks_before_studyentry, current_sact,
         starts_with('signature'),
         tstart, tstop, vte)

# mut sig stats #####
# read in mutational signatures stats
sig_stats <- fread('./results/mutational_signatures_stats.txt')
  # remove any signatures with < 10 VTE cases to prevent issues with model convergence due to sparse data
sigs <- sig_stats %>% filter(cancer == 'PANCANCER' &
                               vte_cases >= 10) %>% pull(signature) %>% unique

# minimally adjusted cox model ####
results <- list()
for (i in sigs) {
  formula <- as.formula(paste('Surv(tstart, tstop, vte) ~', i, '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
  results[[i]] <- coxph(formula, data = sdata) %>% 
  # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
  tidy(exp=T, conf.int = T) %>% 
  # round estimates to 4 decimals
  mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
  # the output provides the HR for VTE associated with each of the exposure variables in the formula,
  # after adjusting for all the other variables
  # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
  rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
  # filter just for the result relating to the mutational signature
  filter(term == i) %>%  
    # select relevant columns for output
  select(term, HR, l95ci, u95ci, pval) %>%
  # add a column which states the model run
  mutate(model_covars = paste(i, '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
  }

results_df <- do.call(rbind, results) %>% arrange(pval) %>% mutate(fdr_p = p.adjust(pval, method = 'fdr'))
  
write.table(results_df, './results/cox_signatures_minadj.tsv', row.names=F, quote =F, sep = '\t')
  
# Fully adjusted model #####
results <- list()
for (i in sigs) {
  formula <- as.formula(paste('Surv(tstart, tstop, vte) ~', i, 
                              '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + disease_type + 
                              stage_grouped + current_sact + sact_over6weeks_before_studyentry'))
  results[[i]] <- coxph(formula, data = sdata) %>% 
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    filter(term == i) %>%
    rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, HR, l95ci, u95ci, pval)
    }

results_df <- do.call(rbind, results) %>%
  # add a column which states the model run
  mutate(model_covars = 'signature + age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + disease_type + stage_grouped + current_sact + sact_over6weeks_before_studyentry') %>%
  arrange(pval)
  
write.table(results_df, './results/cox_signatures_fulladj.tsv', row.names=F, quote =F, sep = '\t')

# Tumour stratified analysis #####
# adjust for age + sex + PC1-4 + stage_grouped + prev_sact + current_sact
# just look at the signatures with FDR-P < 0.1 in the fully adjusted pancancer analysis
## sig subset ####
sig_subset <- fread('./results/cox_signatures_fulladj.tsv') %>% filter(fdr_p <= 0.1) %>% pull(term) %>% unique

# first create vector of disease_types
tumour_types <- ac %>% 
  # remove 'OTHER' because insufficient VTE cases to run stratified analysis 
  filter(disease_type != 'OTHER') %>%
  # create vector of remaining disease types
  pull(disease_type) %>% unique

results <- list()
# for each disease type run the analysis
for (t in tumour_types) {
  temp_data <- sdata %>% filter(disease_type == t)
  results[[t]] <- list()
  for (i in sig_subset) {
    formula <- as.formula(paste('Surv(tstart, tstop, vte) ~', i, 
                                '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + 
                              stage_grouped + current_sact + sact_over6weeks_before_studyentry'))
    results[[t]][[i]] <- coxph(formula, data = temp_data) %>% 
      # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
      tidy(exp=T, conf.int = T) %>% 
      # round estimates to 4 decimals
      mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
      # the output provides the HR for VTE associated with each of the exposure variables in the formula,
      # after adjusting for all the other variables
      # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
      rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
      # select relevant columns for output
      select(term, HR, l95ci, u95ci, pval) %>% 
      # filter only for the results relating to signatures
      filter(term == i) %>%     # add a column which states the type of cancer
      mutate(cancer = t)
               }}
  
results_df <- lapply(results, function(x) do.call(rbind, x)) %>%
  do.call(rbind,.) %>% 
  # add a column which states the model run
  mutate(model_covars = 'signature + age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + stage_grouped + current_sact + sact_over6weeks_before_studyentry') %>%
  # a few of the cancers have uci95 = Inf as the model did not converge due to sparse data
  # for these cancers the analysis is invalid therefore substitute with NA vals
  mutate(HR = case_when(u95ci == 'Inf' ~ NA, .default = HR),
         l95ci = case_when(u95ci == 'Inf' ~ NA, .default = l95ci),
         pval = case_when(u95ci == 'Inf' ~ NA, .default = pval),
         u95ci = case_when(u95ci == 'Inf' ~ NA, .default = u95ci)
         ) %>%
  arrange(term, pval)

# save
write.table(results_df, './results/cox_signatures_tumour_stratified.tsv', row.names=F, quote =F, sep = '\t')

# Tumour interactions #### 
##  ANOVA #####
results_anova <- list()
for (i in sig_subset) {
  # additive model
  additive_model <- coxph(as.formula(paste('Surv(tstart, tstop, vte) ~', i, '+ disease_type +
                            age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                            stage_grouped + 
                            sact_over6weeks_before_studyentry + current_sact')), 
                          data = sdata)
  # interaction cox model
  interaction_model <- coxph(as.formula(paste0(
    'Surv(tstart, tstop, vte) ~',
    i, ' + disease_type:', i, 
    ' + disease_type + age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + stage_grouped + sact_over6weeks_before_studyentry + current_sact')), 
                             data = sdata) 
  results_anova[[i]] <- anova(additive_model, interaction_model)  
}
# extract anova pvalue resutls to table and save
results_anova_df <- tibble(
  signature = names(results_anova),
  p_value = sapply(results_anova, function(x) x$`P(>|Chi|)`[2])
) %>%
  mutate(test = paste('anova(additive model with disease_type +', 
                      signature, 'vs interaction model incorporating disease_type :', signature))

write.table(results_anova_df, './results/anova_tumour_signature_interactions.tsv',
            row.names = F, quote = F, sep = '\t')
rm(results_anova_df, results_anova, interaction_model, additive_model)


# Time interactions ####
# create vector of unique timepoints to split over:
cut1 <- 365*1
cut2 <- 365*2
cut3 <- 365*3
cut4 <- 365*4
# see Therneau et al for other methods


# Now run the model with this time interaction term
res <- list()
for (i in sig_subset) {
  # split each record according to timecuts
  formula <- as.formula(paste('Surv(time_to_status, binary_status) ~', i, 
                              '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
  newdata <- survSplit(formula, data=tdata, 
                       cut=c(cut1, cut2, cut3, cut4), id = 'participant_id') %>%
    # create a new variable time_cut
    mutate(time_point = as.factor(case_when(time_to_status <= cut1 ~ 1,
                                            time_to_status <= cut2 ~ 2,
                                            time_to_status <= cut3 ~ 3,
                                            time_to_status <= cut4 ~ 4,
                                            time_to_status > cut4 ~ 5)))
  
  # now run the model with the linear time interaction term
  new_formula <- as.formula(paste0('Surv(tstart, time_to_status, binary_status) ~ ', i, 
                                   ' + time_point + time_point:', i, '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
  res[[i]] <- coxph(new_formula, data = newdata) %>%
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for sig (after adjusting for age, sex and pcs 1-4 and disease_type)
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    dplyr::rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    dplyr::select(term, HR, l95ci, u95ci, pval)
}

res_df <- do.call(rbind, res) %>% 
    filter(str_detect(term, str_c(sig_subset, collapse = '|'))) %>%
  # be explicit about which model was run
  mutate(model = 'signature + time_point + time_point:signature + age_study_entry + sex + pc1 + pc2 + pc3 + pc4')

write.table(res_df, './results/cox_signatures_time_interaction.tsv', row.names = F, quote =F, sep = '\t')

# time stratified ####
res <- list()
for (i in sig_subset) {
  # split each record according to timecuts
  formula <- as.formula(paste('Surv(time_to_status, binary_status) ~', i, 
                              '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
  newdata <- survSplit(formula, data=tdata, 
                       cut=c(cut1, cut2, cut3, cut4), id = 'participant_id') %>%
    # create a new variable time_cut
    mutate(time_point = as.factor(case_when(time_to_status <= cut1 ~ 1,
                                            time_to_status <= cut2 ~ 2,
                                            time_to_status <= cut3 ~ 3,
                                            time_to_status <= cut4 ~ 4,
                                            time_to_status > cut4 ~ 5)))
  for (t in 1:5) { 
    newdata_timepoint <- newdata %>% filter(time_point == t)
    
    # run the model with the time interaction term, with the data stratified by time point
    new_formula <- as.formula(paste0('Surv(tstart, time_to_status, binary_status) ~ ', i, 
                                     '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
    
    res[[i]][[t]] <- coxph(new_formula, 
                           data=newdata_timepoint) %>%
      tidy(exp=T, conf.int = T) %>% 
      # round estimates to 4 decimals
      mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
      # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
      dplyr::rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
      # select relevant columns for output
      dplyr::select(term, HR, l95ci, u95ci, pval) %>%
      # add a column which states the model run
      mutate(model_covars = paste(i, '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4')) %>%
      # add the timepoint being queried
      mutate(time_point = t) %>%
      # add number of VTE cases occurring in that time period
      mutate(ncase = newdata_timepoint %>% filter(binary_status == 1) %>% nrow)
  }
}

# merge into single df
res_df <- lapply(res, function(x) do.call(rbind, x)) %>%
  do.call(rbind,.) %>%
  filter(str_detect(term, str_c(sig_subset, collapse = '|'))) %>%
  # arrange in order of term and then timepoint
  arrange(term, time_point) %>%
  #reorder columns
  relocate(term, time_point, ncase)

#save
write.table(res_df, './results/cox_signatures_time_stratified.tsv', row.names = F, quote =F, sep = '\t')


# Ancestry-stratification #####
ancestries <- ac %>% 
  # create vector of genetically inferred ancestries (from PCs)
  pull(genetically_inferred_ancestry_nothr) %>% unique

results <- list()
for (a in ancestries) {
  # filter the analysis cohort
  temp_data <- tdata %>% filter(genetically_inferred_ancestry_nothr == a)
  results[[a]] <- list()
  for (i in sig_subset) {
    formula <- as.formula(paste('Surv(time = time_to_status, event = binary_status) ~', i, 
                                '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
     results[[a]][[i]] <- coxph(formula,
                               data = temp_data) %>%
      # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
      tidy(exp=T, conf.int = T) %>% 
      # round estimates to 4 decimals
      mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
      # the output provides the HR for VTE associated with each of the exposure variables in the formula,
      # after adjusting for all the other variables
      # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4 and disease_type)
      filter(str_detect(term, 'signature')) %>% 
      # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
      rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
      # select relevant columns for output
      select(term, HR, l95ci, u95ci, pval) %>%
      mutate(ancestry = a)
  }}

results_df <- lapply(results, function(x) do.call(rbind, x)) %>%
  do.call(rbind,.) %>%
  arrange(ancestry) %>%
  # add a column which states the model run
  mutate(model_covars = 'signature + age_study_entry + sex + pc1-4') %>%
  # a few of the subgroups have uci95 = Inf as the model did not converge due to sparse data
  # for these rows the analysis is invalid therefore substitute with NA vals
  mutate(HR = case_when(u95ci == 'Inf' ~ NA, .default = HR),
         l95ci = case_when(u95ci == 'Inf' ~ NA, .default = l95ci),
         pval = case_when(u95ci == 'Inf' ~ NA, .default = pval),
         u95ci = case_when(u95ci == 'Inf' ~ NA, .default = u95ci)
  )


write.table(results_df, './results/cox_signatures_ancestry_stratified.tsv', row.names=F, quote = F, sep = '\t')
rm(results, results_df, temp_data, a, i)

# Anti-coagulation-stratification #####

results <- list()

for (a in c('yes', 'no')) {
  temp_data  <- tdata %>% filter(prior_anticoag_indication == a)
  results[[a]] <- list()
  for (i in sig_subset) {
    formula <- as.formula(paste('Surv(time = time_to_status, event = binary_status) ~', i, 
                                '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
    results[[a]][[i]] <- coxph(formula,
                               data = temp_data) %>%
      # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
      tidy(exp=T, conf.int = T) %>% 
      # round estimates to 4 decimals
      mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
      # the output provides the HR for VTE associated with each of the exposure variables in the formula,
      # after adjusting for all the other variables
      # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4 and disease_type)
      filter(str_detect(term, 'signature')) %>% 
      # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
      rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
      # select relevant columns for output
      select(term, HR, l95ci, u95ci, pval) %>%
      mutate(prior_anticoag_indication = a)
  }}

# note had issues with model convergence for some tumour-gene combinations 
results_df <- lapply(results, function(x) do.call(rbind, x)) %>%
  do.call(rbind,.) %>%
  # add a column which states the model run
  mutate(model_covars = 'signature + age_study_entry + sex + pc1-4')

write.table(results_df, './results/cox_signatures_anticoag_stratified.tsv', row.names=F, quote = F, sep = '\t')
rm(results, results_df, temp_data, a, i)

# New-untreated cancer -stratification #####
ac_sensitivity_diag_equals_time0 <- fread('./data/ac_sensitivity_diag_equals_time0.csv') %>%
  # create a binary status column for cox models where death is censored
  mutate(binary_status = case_when(status == 1 ~ 1, .default = 0),
            age_study_entry = as.integer(format(study_entry, "%Y")) - year_of_birth)
# get sigs for this c ohort
table_queried <- "cancer_analysis" 
signatures <- paste0("signature_", 1:30)
columns_selected <- c("participant_id", "tumour_sample_platekey", signatures)

tmb_everyone <- labkey_select() %>%
  # convert signatures (continuous numerical with high proprtion NA/zero vals to a binary val)
  # if NA then code as zero, otherwise code as 1
   mutate(across(starts_with('signature'), ~ifelse(is.na(.), 0,1)))

rm(table_queried, signatures, columns_selected)
# subset the new untreated cancers from the ac  
df_new_untreated <- ac_sensitivity_diag_equals_time0 %>% 
  filter(new_untreated_cancer == 'yes') %>%
  # add sig info
  left_join(., tmb_everyone, by = join_by(participant_id, tumour_sample_platekey))
rm(ac_sensitivity_diag_equals_time0)
# check they all have signatures
table(df_new_untreated$signature_8)

results <- list()
for (i in sig_subset) {
    formula <- as.formula(paste('Surv(time = time_to_status, event = binary_status) ~', i, 
                                '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
        results[[i]] <- coxph(formula,
                               data = df_new_untreated) %>%
      # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
      tidy(exp=T, conf.int = T) %>% 
      # round estimates to 4 decimals
      mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
      # the output provides the HR for VTE associated with each of the exposure variables in the formula,
      # after adjusting for all the other variables
      filter(str_detect(term, 'signature')) %>% 
      # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
      rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
      # select relevant columns for output
      select(term, HR, l95ci, u95ci, pval) %>%
    mutate(model_covars = paste(term, '+ age_study_entry + sex + pc1-4'),
    sensitivity_analysis = 'new-untreated; time zero as cancer diagnosis')
  }

# note had issues with model convergence for some tumour-gene combinations 
results_df <- do.call(rbind,results)

write.table(results_df, './results/cox_new_untreated_timezero_diag_signatures.tsv', row.names=F, quote = F, sep = '\t')

# immune check point inhibition #####

# load in ac_ICI_CDKi
ac_ICI_CDKi <- fread('./data/ac_ICI_CDKi.csv')

# total no participants receiving ICI
ac_ICI_CDKi %>% filter(!is.na(time_to_ICI)) %>% count

# create a new tdata frame which contains the ICI info
tdata_ICI <- left_join(tdata, ac_ICI_CDKi, by = 'participant_id') %>% 
  # if ICI given within -42 days of study entry set time as 0
  mutate(time_to_ICI = case_when(time_to_ICI <= 0 &
                                   time_to_ICI >= -42 ~ 0,
                                 # if ICI given more than -42 days before study entry set time as NA                          
                                 time_to_ICI < -42 ~ NA,
                                 TRUE ~ time_to_ICI))


sdata_ICI <- tmerge(data1 = tdata_ICI, data2 = tdata_ICI,
                    vte = event(time_to_status, binary_status), 
                    current_sact = tdc(time_to_current_sact),
                    ICI = tdc(time_to_ICI),
                    id = participant_id,
                    options= list(idname="participant_id")) %>%
  # select only relevant columns
  select(participant_id, tumour_sample_platekey, age_study_entry, age_study_entry_scaled, age_cancerdiag, 
         sex, pc1, pc2, pc3, pc4,
         disease_type, stage_numeric_imputed, stage_grouped,
         sact_over6weeks_before_studyentry, current_sact, time_to_ICI,ICI,
         starts_with('signature'),
         tstart, tstop, vte)

# now run different models for each signature with FDR_P < 0.1.
sig_subset <- fread('./results/cox_signatures_fulladj.tsv') %>% filter(fdr_p < 0.1) %>% pull(term) %>% unique

results <- list()
for (s in sig_subset) {
  results[[s]] <- list()# create a sublist
# create the models
full_adj <- 'sact_over6weeks_before_studyentry + current_sact'
#2) interaction between sigature:current_sact
sact_interaction <- paste0('current_sact + ', s, ':current_sact')
# 2) the model which contains just ICI
ICI <- 'ICI'
# 3) model which looks at interaction between ICI and signature
ICI_interaction <- paste0('ICI + ',s,':ICI')
# create a list of all the models
models <- c(full_adj, sact_interaction, ICI, ICI_interaction)
 # now run each model
  for (f in models) {
  formula <- as.formula(paste('Surv(tstart, tstop, vte) ~', s,'+ 
                                age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                                disease_type + stage_grouped +', f))
  results[[s]][[f]] <- coxph(formula, data = sdata_ICI) %>% 
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    filter(str_detect(term, 'signature')|str_detect(term, 'ICI')|str_detect(term, 'sact')) %>%
    # select relevant columns for output
    select(term, HR, l95ci, u95ci, pval) %>%
    mutate(model = paste(s,'+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + disease_type + stage_grouped +', f))
}}

results_df <- lapply(results, function(x) do.call(rbind, x)) %>% do.call(rbind, .)

write.table(results_df, './results/cox_signature_ICI_interactions.tsv', row.names=F, quote =F, sep = '\t')

## Continuous model ####
# model just the continuous non-zero values
tdata_continuous <- ac %>% 
  # merge with the TMB data with the CONTINUOUS sig values
  left_join(., tmb_ac_continuous, by = 'tumour_sample_platekey') %>%
  # add +0.5 to time_to_status only where time_to_status =0
  # this step is required becuse tmerge does not accept stop times of 0
  # see vignette by Therneau et al
  mutate(time_to_status = case_when(time_to_status == 0 ~ 0.5,
                                    .default = time_to_status))
results <- list()
for (i in sig_subset) {
  # subset sdata_sensitivity to exclude outliers
  temp_data <- tdata_continuous %>% filter(!is.na(!!sym(i)))
  formula <- as.formula(paste('Surv(time_to_status, binary_status) ~', i, 
                              '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
  results[[i]] <- coxph(formula, data = temp_data) %>% 
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    filter(str_detect(term, i)) %>%
    rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, HR, l95ci, u95ci, pval) %>%
    mutate(model_covars = paste(i, '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + disease_type + stage_grouped + current_sact + sact_over6weeks_before_studyentry'))
  
}

results_df <- do.call(rbind, results)
write.table(results_df, './results/cox_signature_continuous_nonzeros.tsv', row.names=F, quote =F, sep = '\t')

# End ####