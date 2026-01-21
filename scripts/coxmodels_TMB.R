# Cox regression models for TMB analyses
# N Cornish

# Config ####
# load config script
setwd(Sys.getenv("wd"))
source('./scripts/config421.R')

# Load formatted analysis cohort ####
ac <- readRDS('./data/analysis_cohort_formatted.rds')

# load the file with the TMBs and mutational signatures 
tmb_ac <- fread('./data/tmb_and_mutsig_ac.csv') %>% 
  # select only the columns relating to TMB
  select(tumour_sample_platekey, somatic_coding_variants_per_mb, starts_with('tmb')) %>%
   # create a column tmb_binary which is high if somatic_coding_mutations_per_mb >= 20
  mutate(tmb_binary_thresh20 = factor(case_when(
    somatic_coding_variants_per_mb < 20 ~ 'tmb_below_20',
    somatic_coding_variants_per_mb >= 20 ~ 'tmb_above_20'
  ), levels = c('tmb_below_20', 'tmb_above_20'))) %>%
  # also create a column tmb_categories which breaks TMB into blocks of 5
  mutate(tmb_categories = factor(case_when(
    somatic_coding_variants_per_mb < 5 ~ 'tmb_below_5',
    somatic_coding_variants_per_mb >= 5 & somatic_coding_variants_per_mb < 10 ~ 'tmb_between_5_10',
    somatic_coding_variants_per_mb >= 10 & somatic_coding_variants_per_mb < 20 ~ 'tmb_between_10_20',
    somatic_coding_variants_per_mb >= 20  ~ 'tmb_above_20',
                                       .default = NA), 
    levels = c('tmb_below_5','tmb_between_5_10','tmb_between_10_20', 'tmb_above_20')))
                                         
# create sdata where sact is coded as a time varying covariate
##tmerge 
# use tmerge because sact is a time dependent covar
# first create the basic dataset for tmerge
tdata <- ac %>% 
  # merge with the TMB data
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
         somatic_coding_variants_per_mb, starts_with('tmb'),
         tstart, tstop, vte)


# run the cox models
## tmb_variables ####
tmb_variables <- sdata %>% 
  select(starts_with('tmb')) %>% names %>% unlist


# minimally adjusted cox model ####
results <- list()
for (i in tmb_variables) {
  formula <- as.formula(paste('Surv(tstart, tstop, vte) ~', i, '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
  results[[i]] <- coxph(formula, data = sdata) %>% 
  # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
  tidy(exp=T, conf.int = T) %>% 
  # round estimates to 4 decimals
  mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
  # the output provides the HR for VTE associated with each of the exposure variables in the formula,
  # after adjusting for all the other variables
  # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
  dplyr::rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
  # select relevant columns for output
  select(term, HR, l95ci, u95ci, pval) %>%
  # add a column which states the model run
  mutate(model_covars = paste(i,'+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
  }

results_df <- do.call(rbind, results) %>%
  # select the rows which relate to TMB not the other covars
  filter(str_detect(term, 'tmb'))
  
write.table(results_df, './results/cox_TMB_minadj.tsv', row.names=F, quote =F, sep = '\t')
  
# Fully adjusted model #####
results <- list()
for (i in tmb_variables) {
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
    dplyr::rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, HR, l95ci, u95ci, pval)
    }

results_df <- do.call(rbind, results) %>%
    # select the rows which relate to TMB not the other covars
  filter(str_detect(term, 'tmb'))   
write.table(results_df, './results/cox_TMB_fulladj.tsv', row.names=F, quote =F, sep = '\t')

## martingale residuals for raw TMB  ####
i <- 'somatic_coding_variants_per_mb'
formula <- as.formula(paste('Surv(tstart, tstop, vte) ~', i, 
                            '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + disease_type + 
                              stage_grouped + current_sact + sact_over6weeks_before_studyentry'))
model <- coxph(formula, data = sdata)
mart_resid <- residuals(model, type = "martingale")
# Create plotting data frame
plot_data <- data.frame(
  TMB = sdata$somatic_coding_variants_per_mb,
  residuals = mart_resid,
  status_binary = factor(sdata$vte, 
                         levels = c(0, 1),
                         labels = c("Censored", "VTE")))

plot <- ggplot(plot_data, aes(x = TMB, 
                      y = residuals, color = status_binary)) +
  geom_point() +
  geom_smooth(method = "loess", span = 2/3, se = TRUE, color = "black", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Censored" = "blue", "VTE" = "red")) +
  labs(x = "TMB (somatic_coding_variants_per_mb)", 
       y = "Martingale Residuals",
       title = "Martingale Residuals vs TMB",
       color = "VTE Status") +
  theme_minimal() +
  theme(legend.position = "bottom")

png('./results/TMB_martingale_residuals.png', width = 800, height= 600)
plot
dev.off()
# data not shown

# Tumour stratified analysis #####
# adjust for age + sex + PC1-4 + stage_grouped + prev_sact + current_sact

## tumour_types ####
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
  for (i in tmb_variables) {
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
      # filter only for the results relating to tmb
      filter(str_detect(term, 'tmb')) %>%
             # add a column which states the type of cancer
      mutate(cancer = t)
               }}
  
results_df <- lapply(results, function(x) do.call(rbind, x)) %>%
  do.call(rbind,.) %>% 
  # a few of the cancers have uci95 = Inf as the model did not converge for tmb_binary due to sparse data
  # for these cancers the analysis is invalid therefore substitute with NA vals
  mutate(HR = case_when(u95ci == 'Inf' ~ NA, .default = HR),
         l95ci = case_when(u95ci == 'Inf' ~ NA, .default = l95ci),
         pval = case_when(u95ci == 'Inf' ~ NA, .default = pval),
         u95ci = case_when(u95ci == 'Inf' ~ NA, .default = u95ci)
  ) %>%
  arrange(term, pval) %>%
  # add a column which states the model run
  mutate(model_covars = 'tmb + age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + stage_grouped + current_sact + sact_over6weeks_before_studyentry')


# save
write.table(results_df, './results/cox_TMB_tumour_stratified.tsv', row.names=F, quote =F, sep = '\t')

# Tumour interactions #### 
##  ANOVA #####
# additive model
additive_model <- coxph(Surv(tstart, tstop, vte) ~ tmb_binary_thresh20 + disease_type +
                            age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                            stage_grouped + 
                            sact_over6weeks_before_studyentry + current_sact, 
                          data = sdata)
# interaction cox model
interaction_model <- coxph(Surv(tstart, tstop, vte) ~ tmb_binary_thresh20
                           + disease_type + disease_type:tmb_binary_thresh20 +
                               age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                              stage_grouped + 
                               sact_over6weeks_before_studyentry + current_sact, 
                             data = sdata) 
results_anova <- anova(additive_model, interaction_model)  

# extract anova pvalue resutls to table and save
results_anova_df <- tibble(
  test = 'anova comparing TMB>=20:disease_type interaction model with additive model',
  p_value = results_anova$`P(>|Chi|)`[2])

write.table(results_anova_df, './results/anova_tumour_TMB_interactions.tsv',
            row.names = F, quote = F, sep = '\t')
rm(results_anova_df, results_anova, interaction_model, additive_model)


# TMB summary stats ####
#look at the range of TMB for each tumour
# use tdata to avoid double-counting patients who are split across multiple rows in sdata
results <- list()
for (t in tumour_types) {
  tmb <- tdata %>% 
    # filter for tumour type
    filter(disease_type == t) %>%
    # ensure each participant represented only once
    distinct(tumour_sample_platekey, .keep_all = T) %>%
    # look at the TMB ranges for that tumour
    select(somatic_coding_variants_per_mb) %>% pull
  high_tmb <- tmb[tmb >= 10] %>% length
  very_high_tmb <- tmb[tmb >= 20] %>% length
  # summary stats
  results[[t]] <- data.frame(tumour = t,
                             tmb_min = min(tmb),
                             q1 =  quantile(tmb, probs = c(0.25)),
                             tmb_median = median(tmb),
                             q3 =  quantile(tmb, probs = c(0.75)),
                             tmb_max = max(tmb),
                             tmb_high_prop = round(high_tmb/length(tmb), digits =3),
                             tmb_very_high_prop = round(very_high_tmb/length(tmb), digits =3)) }

# other was not in the tumour_types vector so run code separately to include
tmb <- tdata %>% filter(disease_type == 'OTHER') %>% distinct(tumour_sample_platekey, .keep_all = T) %>% select(somatic_coding_variants_per_mb) %>% pull
high_tmb <- tmb[tmb >= 10] %>% length
very_high_tmb <- tmb[tmb >= 20] %>% length
results[['OTHER']] <- data.frame(tumour = 'OTHER',
                                 tmb_min = min(tmb),
                                 q1 =  quantile(tmb, probs = c(0.25)),
                                 tmb_median = median(tmb),
                                 q3 =  quantile(tmb, probs = c(0.75)),
                                 tmb_max = max(tmb),
                                 tmb_high_prop = round(high_tmb/length(tmb), digits =3),
                                 tmb_very_high_prop = round(very_high_tmb/length(tmb), digits =3)) 
# PANCANCER result
tmb <- tdata %>% distinct(tumour_sample_platekey, .keep_all = T) %>% select(somatic_coding_variants_per_mb) %>% pull
high_tmb <- tmb[tmb >= 10] %>% length
very_high_tmb <- tmb[tmb >= 20] %>% length
results[['PANCANCER']] <- data.frame(tumour = 'PANCANCER',
                                     tmb_min = min(tmb),
                                     q1 =  quantile(tmb, probs = c(0.25)),
                                     tmb_median = median(tmb),
                                     q3 =  quantile(tmb, probs = c(0.75)),
                                     tmb_max = max(tmb),
                                   tmb_high_prop = round(high_tmb/length(tmb), digits =3),
                                   tmb_very_high_prop = round(very_high_tmb/length(tmb), digits =3)) 
# combine into single df
results_df <- do.call(rbind, results) %>% # round to 2 decimals
  mutate_if(is.numeric, ~round(., digits = 3))

write.table(results_df, './results/TMB_summarystats.tsv', row.names=F, quote =F, sep = '\t')


# create a separate boxplot of TMB by cancer type
library(ggplot2)

# create a factor of cancer-typesin order of median TMB
tmb_summstats <- fread('./results/TMB_summarystats.tsv')
tmb_cancer_factor <- tmb_summstats %>% 
  #remove pancancer
  filter(tumour != 'PANCANCER') %>% 
  # arrange in order of TMB
  arrange(desc(tmb_median)) %>%
  pull(tumour) %>% factor(as.character(.), levels=unique(.))


# boxplot
png('./results/TMB_raw_distribution_with_outliers.png', width=800, height=600)
ggplot(data = tdata %>%
         # convert disease type to the order in tmb_cancer_factor
         mutate(disease_type = factor(disease_type, levels = tmb_cancer_factor)) %>%
         # remove any zero values for TMB otherwise can't plot y axis on log2 scale
         filter(somatic_coding_variants_per_mb > 0),
       # plot disease type on x axis and TMB on y axis  
       aes(x = disease_type, y = somatic_coding_variants_per_mb)) +
  # remove outlying points but plot full range of vals
  geom_boxplot()+
  coord_trans(y = "log2") +
  scale_y_continuous(breaks = c(0.12, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512)) +
  labs(title = "TMB Distribution by Cancer Type (with outliers)",
       x = "Cancer Type",
       y = "raw TMB \n (coding mutations per Mb)") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        # create a margin otherwise text spreads off the edge of graph
        plot.margin = margin(10, 10, 10, 80))
dev.off()

# Time interactions ####
# check  proportional hazard assumptions in the MINIMALLY ADJUSTED ANALYSIS
# for this I need to store the complete CoxPH model
# use tdata for simplicity (don't need to split by SACT time)

# here I include gene:time_to_status in the model (linear interaction term between the gene and time)
# create vector of unique timepoints to split over:
cut1 <- 365*1
cut2 <- 365*2
cut3 <- 365*3
cut4 <- 365*4

res <- list()
for (i in tmb_variables) {
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
    # extract just the HR for TMB (after adjusting for age, sex and pcs 1-4 and disease_type)
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    dplyr::rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    dplyr::select(term, HR, l95ci, u95ci, pval)
}

res_df <- do.call(rbind, res) %>% 
  filter(str_detect(term, str_c(tmb_variables, collapse = '|'))) %>%
  # be explicit about which model was run
  mutate(model = 'tmb + time_point + time_point:tmb + age_study_entry + sex + pc1 + pc2 + pc3 + pc4')


write.table(res_df, './results/cox_TMB_time_interaction.tsv', row.names = F, quote =F, sep = '\t')

# time stratified ####
res <- list()
for (i in tmb_variables) {
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
      # the output provides the HR for VTE associated with each of the exposure variables in the formula,
      # after adjusting for all the other variables
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

# merge into a single df
res_df <- lapply(res, function(x) do.call(rbind, x)) %>%
  do.call(rbind,.) %>%
  filter(str_detect(term, str_c(tmb_variables, collapse = '|'))) %>%
  # arrange in order of term and then timepoint
  arrange(term, time_point) %>%
  #reorder columns
  relocate(term, time_point, ncase)

#save
write.table(res_df, './results/cox_TMB_time_stratified.tsv', row.names = F, quote =F, sep = '\t')

# Ancestry-stratification #####
ancestries <- ac %>% 
  # create vector of genetically inferred ancestries (from PCs)
  pull(genetically_inferred_ancestry_nothr) %>% unique

results <- list()
# for each disease type run the analysis
for (a in ancestries) {
  # filter the analysis cohort 
  temp_data <- tdata %>% filter(genetically_inferred_ancestry_nothr == a)
  # create a new list to store results
  results[[a]] <- list()
  for (i in tmb_variables) {
    formula <- as.formula(paste('Surv(time = time_to_status, event = binary_status) ~', i, 
                                '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
    # now run the cox model adjusting for everything except tumour type
    results[[a]][[i]] <- coxph(formula,
                               data = temp_data) %>%
      # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
      tidy(exp=T, conf.int = T) %>% 
      # round estimates to 4 decimals
      mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
      # the output provides the HR for VTE associated with each of the exposure variables in the formula,
      # after adjusting for all the other variables
      filter(str_detect(term, 'somatic')|str_detect(term, 'tmb')) %>% 
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
  mutate(model_covars = 'TMB + age_study_entry + sex + pc1-4') %>%
  # a few of the subgroups have uci95 = Inf as the model did not converge due to sparse data
  # for these rows the analysis is invalid therefore substitute with NA vals
  mutate(HR = case_when(u95ci == 'Inf' ~ NA, .default = HR),
         l95ci = case_when(u95ci == 'Inf' ~ NA, .default = l95ci),
         pval = case_when(u95ci == 'Inf' ~ NA, .default = pval),
         u95ci = case_when(u95ci == 'Inf' ~ NA, .default = u95ci)
  )


write.table(results_df, './results/cox_TMB_ancestry_stratified.tsv', row.names=F, quote = F, sep = '\t')
rm(results, results_df, temp_data, a, i)
# Anti-coagulation-stratification #####

results <- list()
for (a in c('yes', 'no')) {
  # filter the analysis cohort for tumour type t
  temp_data  <- tdata %>% filter(prior_anticoag_indication == a)
  # create a new list to store results for that tumour
  results[[a]] <- list()
   for (i in tmb_variables) {
    formula <- as.formula(paste('Surv(time = time_to_status, event = binary_status) ~', i, 
                                '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4'))
    # now run the cox model adjusting for everything except tumour type
    results[[a]][[i]] <- coxph(formula,
                               data = temp_data) %>%
      # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
      tidy(exp=T, conf.int = T) %>% 
      # round estimates to 4 decimals
      mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
      # the output provides the HR for VTE associated with each of the exposure variables in the formula,
      # after adjusting for all the other variables
      # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4 and disease_type)
      filter(str_detect(term, 'somatic')|str_detect(term, 'tmb')) %>% 
      # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
      rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
      # select relevant columns for output
      select(term, HR, l95ci, u95ci, pval) %>%
      mutate(prior_anticoag_indication = a)
  }}


results_df <- lapply(results, function(x) do.call(rbind, x)) %>%
  do.call(rbind,.) %>%
  # add a column which states the model run
  mutate(model_covars = 'TMB + age_study_entry + sex + pc1-4')

write.table(results_df, './results/cox_TMB_anticoag_stratified.tsv', row.names=F, quote = F, sep = '\t')
rm(results, results_df, temp_data, a, i)

# New-untreated cancer -stratification #####
ac_sensitivity_diag_equals_time0 <- fread('./data/ac_sensitivity_diag_equals_time0.csv') %>%
  # create a binary status column for cox models where death is censored
  mutate(binary_status = case_when(status == 1 ~ 1, .default = 0),
         # calculate age study entry 
         age_study_entry = as.integer(format(study_entry, "%Y")) - year_of_birth)

# pull TMB for new ac
table_queried <- "cancer_analysis"
columns_selected <- c("participant_id", "tumour_sample_platekey",
                      "somatic_coding_variants_per_mb")
tmb_everyone <- labkey_select() %>% 
  # create the binary tmb_threshold
  mutate(tmb_binary_thresh20 = factor(case_when(
    somatic_coding_variants_per_mb < 20 ~ 'tmb_below_20',
    somatic_coding_variants_per_mb >= 20 ~ 'tmb_above_20'
  ), levels = c('tmb_below_20', 'tmb_above_20')))

# subset the new untreated cancers from the ac  
df_new_untreated <- ac_sensitivity_diag_equals_time0 %>% 
  filter(new_untreated_cancer == 'yes') %>%
  # add TMB info
  left_join(., tmb_everyone, by = join_by(participant_id, tumour_sample_platekey))
rm(ac_sensitivity_diag_equals_time0)
 # check they all have TMB
summary(df_new_untreated$tmb_binary_thresh20)


results  <- coxph(Surv(time = time_from_diagnosis_to_status, 
                       event = binary_status) ~ tmb_binary_thresh20 + age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                      data = df_new_untreated) %>% 
  # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
  tidy(exp=T, conf.int = T) %>% 
  # round estimates to 4 decimals
  mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
  # the output provides the HR for VTE associated with each of the exposure variables in the formula,
  # after adjusting for all the other variables
  # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
  filter(str_detect(term, 'tmb_binary_thresh20')) %>% 
  # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
  rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
  # select relevant columns for output
  select(term, HR, l95ci, u95ci, pval) %>%
  # add a column which states the model run
  mutate(model_covars = paste(term, '+ age_study_entry + sex + pc1-4'),
         sensitivity_analysis = 'new-untreated; time zero as cancer diagnosis')


write.table(results, './results/cox_new_untreated_timezero_diag_TMB.tsv', row.names=F, quote = F, sep = '\t')

# immune check point inhibition #####

# load in ac_ICI_CDKi
ac_ICI_CDKi <- fread('./data/ac_ICI_CDKi.csv')

# total no participants receiving ICI
ac_ICI_CDKi %>% filter(!is.na(time_to_ICI)) %>% count

# create a new tdata frame which contains the ICI info
tdata_ICI <- left_join(tdata, ac_ICI_CDKi, by = 'participant_id') %>% 
  # if CDKi given within -42 days of study entry set time as 0
  mutate(time_to_ICI = case_when(time_to_ICI <= 0 &
                                   time_to_ICI >= -42 ~ 0,
                                 # if CDKi given more than -42 days before study entry set time as NA                          
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
         tmb_binary_thresh20,
         tstart, tstop, vte)

# now run different models. 
#1) the normal 'fully adjusted' model with any_sact
full <- c('sact_over6weeks_before_studyentry + current_sact')
#2) interaction between any_sact and tmb
sact_interaction <- c('current_sact + tmb_binary_thresh20:current_sact')
# 2) the model which contains just ICI
ICI <- c('ICI')
# 3) model which looks at interaction between ICI and tmb
ICI_interaction <- c('ICI + tmb_binary_thresh20:ICI')
# all models together
models <- c(full, sact_interaction, ICI, ICI_interaction)

# just run for TMB binary thresh20
results <- list()
for (f in models) {
  formula <- as.formula(paste('Surv(tstart, tstop, vte) ~ tmb_binary_thresh20 + 
                                age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                                disease_type + stage_grouped +', f))
  results[[f]] <- coxph(formula, data = sdata_ICI) %>% 
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    filter(str_detect(term, 'tmb_binary_thresh20')|str_detect(term, 'ICI')|str_detect(term, 'sact')) %>%
    # select relevant columns for output
    select(term, HR, l95ci, u95ci, pval) %>%
    mutate(model = paste('tmb_binary +', f))
}

results_df <- do.call(rbind, results)

write.table(results_df, './results/cox_TMB_ICI_interactions.tsv', row.names=F, quote =F, sep = '\t')
# data not shown in manuscript

# compare the model with vs without the interaction terms
formula <- as.formula(paste('Surv(tstart, tstop, vte) ~ tmb_binary_thresh20 + 
                                age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                                disease_type + stage_grouped + ICI'))
additive_model <- coxph(formula, data = sdata_ICI)

formula2 <- as.formula(paste('Surv(tstart, tstop, vte) ~ tmb_binary_thresh20 + 
                                age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                                disease_type + stage_grouped + ICI + tmb_binary_thresh20:ICI'))
interaction_model <- coxph(formula2, data = sdata_ICI)

# compare the two models
sink('./results/TMB_ICI_anova.txt')
anova(additive_model, interaction_model, test = 'LRT') %>% print
sink()

## leave_one_out_analysis ####
results <- list()
for (i in (ac %>% select(disease_type) %>% pull %>% unique)) {
  sdata_leave_one_out <- sdata %>% filter(disease_type != i)
  formula <- as.formula('Surv(tstart, tstop, vte) ~ tmb_binary_thresh20 + age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + disease_type + 
                              stage_grouped + current_sact + sact_over6weeks_before_studyentry')
  results[[i]] <- coxph(formula, data = sdata_leave_one_out) %>% 
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    mutate(cancer_left_out = i) %>%
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    dplyr::rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, HR, l95ci, u95ci, pval, cancer_left_out)
}
sdata_leave_one_out <- sdata  %>%
  mutate(tmb_binary_over20 = factor(case_when(somatic_coding_variants_per_mb < 20 ~ 'below_20',
                                              somatic_coding_variants_per_mb >= 20 ~ 'above_20'), 
                                    levels = c('below_20', 'above_20')))

results[['PANCANCER']] <- coxph(formula, data = sdata_leave_one_out) %>% 
  # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
  tidy(exp=T, conf.int = T) %>% 
  mutate(cancer_left_out = 'PANCANCER') %>%
  # round estimates to 4 decimals
  mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
  # the output provides the HR for VTE associated with each of the exposure variables in the formula,
  # after adjusting for all the other variables
  # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
  dplyr::rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
  # select relevant columns for output
  select(term, HR, l95ci, u95ci, pval, cancer_left_out)
results_df <- do.call(rbind, results) %>%
    # select the rows which relate to TMB not the other covars
  filter(str_detect(term, 'tmb_binary'))   

write.table(results_df, './results/TMB_leave_one_out.tsv', quote = F, row.names =F)


# End ####