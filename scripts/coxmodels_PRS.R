# Cox models PRS (polygenic risk score)
# script for evaluating effect of germline PRS on VTE risk
# N Cornish

setwd(Sys.getenv('wd'))
source('./scripts/config421.R')
library(survival)
# Load formatted analysis cohort ####
ac <- readRDS('./data/analysis_cohort_formatted.rds')

#Load in the germline polygenic risk scores from Klarin et al for all participants 
# (extended 297(+2) SNP PRS = 299 SNPs total). 
# See README and plink_calculate_PRS.sh for details on how this was calculated
prs_scores <- fread('./data/participant_PRS_scores_klarin.txt.sscore') %>%
  # filter only for the IID column and SCORE1_AVG column. Note the average means the total PRS was then divided by the number of valid allele observations for that participant - this denominator was 299 if a participant had a valid genotype for every single SNP but due to quality control most participants had genotypes missing for a number of the SNPs.
  dplyr::select(IID, SCORE1_AVG) %>%
  # change column names so they are more intuitive
  dplyr::rename('germline_sample_platekey' = IID,
                'klarin_PRS' = SCORE1_AVG) %>%
  # this table contains PRS scores for all participants in AggV2: filter specifically for the analysis cohort
  filter(germline_sample_platekey %in% ac$germline_sample_platekey) %>%
  # scale the PRS to have a mean of zero and SD of 1
  mutate(klarin_PRS_scaled = scale(klarin_PRS) %>% as.vector) %>%
  # create a PRS quartile column
  mutate(prs_quartile = ntile(klarin_PRS_scaled, 4)) %>%
  mutate(prs_category = case_when(prs_quartile == 1 ~ 'lower_3_quartiles',
                                  
                                  prs_quartile == 2 ~ 'lower_3_quartiles',
                                  
                                  prs_quartile == 3 ~ 'lower_3_quartiles',
                                  
                                  prs_quartile == 4 ~ 'top_quartile',
                                  
                                  .default = NA
                                  
  ))

# check 
mean(prs_scores$klarin_PRS_scaled) %>% round
sd(prs_scores$klarin_PRS_scaled)

#count any patients who are missing PRS scores because they are not present in AggV2
ac %>% filter(!germline_sample_platekey %in% prs_scores$germline_sample_platekey) %>% count
# these patients will be excluded from the PRS analysis

# now run the models
# create sdata where sact is coded as a time varying covariate
##tmerge 
tdata <- ac %>% 
  # merge with the TMB data
  left_join(., prs_scores, by = 'germline_sample_platekey') %>%
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
         starts_with('klarin'), starts_with('prs'),
         tstart, tstop, vte)


# run the cox models
# look at associations when PRS categorised by quartile and also when PRS modelled continuously (per SD)
prs_variables <- sdata %>% select(klarin_PRS_scaled, starts_with('prs')) %>% names %>% unlist


# minimally adjusted cox model ####
results <- list()
for (i in prs_variables) {
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
    # select relevant columns for output
    select(term, HR, l95ci, u95ci, pval) %>%
    # add a column which states the model run
    mutate(model_covars = 'prs + age_study_entry + sex + pc1 + pc2 + pc3 + pc4')
}

results_df <- do.call(rbind, results) %>%
  # select the rows which relate to PRS not the other covars
  filter(term == 'klarin_PRS_scaled'|term == 'prs_quartile'|term == 'prs_categorytop_quartile')

write.table(results_df, './results/cox_PRS_minadj.tsv', row.names=F, quote =F, sep = '\t')

# Fully adjusted model #####
results <- list()
for (i in prs_variables) {
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
    rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, HR, l95ci, u95ci, pval)
}

results_df <- do.call(rbind, results) %>%
  mutate(model_covars = 'PRS + age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + disease_type + stage_grouped + current_sact + sact_over6weeks_before_studyentry') %>%
  # select the rows which relate to PRS not the other covars
  filter(term == 'klarin_PRS_scaled'|term == 'prs_quartile'|term == 'prs_categorytop_quartile')

write.table(results_df, './results/cox_PRS_fulladj.tsv', row.names=F, quote =F, sep = '\t')


# Gene interactions ######
# Now look at interactions between the PRS and somatic mutations
all_variants <-fread('./data/combined_variants.csv') %>% 
  # participants will only be counted once for each gene even if they carry multiple different variant
  distinct
fdr_sig_genes <- fread('./results/cox_genes_minadj.tsv') %>% filter(fdr_p <= 0.1) %>% pull(gene)

results <- list()
for (i in fdr_sig_genes) {
  # filter the all_variants table for the gene in question
  variants_by_gene <- all_variants %>% filter(gene == i)
  # create a new column in the analysis cohort dataframe 
  temp_data <- tdata %>% 
    # which identifies whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  # now run the cox model
  results[[i]] <- coxph(Surv(time = time_to_status, event = binary_status) ~ gene + klarin_PRS_scaled + gene:klarin_PRS_scaled +
                          age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                        data = temp_data) %>% 
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
    filter(str_detect(term, 'gene')|str_detect(term, 'klarin')) %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, HR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(term = gsub('gene', i, term)) %>%
    mutate(model_covars = paste0(i, ' + klarin_PRS_scaled + ', i,
                                 ':klarin_PRS_scaled + age_study_entry + sex + pc1-4'))
}

# merge into a single df and look at the fdr adjusted p values
results_df <- do.call(rbind, results) 

write.table(results_df, './results/cox_PRS_gene_interactions.tsv', row.names=F, quote = F, sep = '\t')

## PCDH15:PRS
# Explore interaction between PRS and PCDH15 in more detail
# PCDH15 alone
i <- 'PCDH15'
variants_by_gene <- all_variants %>% filter(gene == i)
# create a new column in the analysis cohort dataframe 
temp_data <- tdata %>% 
  # which identifies whether the participants sample appears in the variants_by_gene table
  mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
# now run the cox model
results_pcdh15 <- coxph(Surv(time = time_to_status, event = binary_status) ~ gene +
                        age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                      data = temp_data) 

results_pcdh15_prs <- coxph(Surv(time = time_to_status, event = binary_status) ~ gene +
                          + klarin_PRS_scaled + age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                        data = temp_data) 

results_pcdh15_prs_int <- coxph(Surv(time = time_to_status, event = binary_status) ~ gene +
                              + klarin_PRS_scaled + gene:klarin_PRS_scaled + age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                            data = temp_data) 
# run anova test to compare the additive model and interaction model
sink('./results/PRS_PCDH15_anova.txt')
anova(results_pcdh15_prs, results_pcdh15_prs_int) %>% print
sink()

# Stratified analyses for PCDH15 based on PRS
temp_data_high_PRS <- temp_data %>% filter(prs_category == 'top_quartile')
high_prs <- coxph(Surv(time = time_to_status, event = binary_status) ~ gene +
                    + age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                  data = temp_data_high_PRS)
temp_data_low_PRS <- temp_data %>% filter(prs_category == 'lower_3_quartiles')
low_prs <- coxph(Surv(time = time_to_status, event = binary_status) ~ gene +
                    + age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                  data = temp_data_low_PRS)

## PCDH15:PRS cuminc curves #####
# split the cohort into high PRS + PCDH15_mutated
# high PRS, No PCDH15 mutation
# low PRS; PCDH15 mutated
# low PRS; no PCDH15 mutation
variants_by_gene <- all_variants %>% filter(gene == 'PCDH15')
# create a new column in the analysis cohort dataframe 
temp_data <- tdata %>% 
  # which identifies whether the participants sample appears in the variants_by_gene table
  mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))

fit_data <- temp_data %>% 
  # split into the 4 categories described above
  mutate(combined_category = 
           case_when(prs_category == 'top_quartile' & gene == 1 ~ 'high PRS; PCDH15 mutated',
                     prs_category == 'top_quartile' & gene == 0 ~ 'high PRS; PCDH15 wiltype',
                     prs_category == 'lower_3_quartiles' & gene == 1 ~ 'low PRS; PCDH15 mutated',
                     prs_category == 'lower_3_quartiles' & gene == 0 ~ 'low PRS; PCDH15 wildtype',
                     .default = NA)
  )

# run a fine-gray weighted model to avoid over-estimating the incidence of VTE
# take 'fit_data' produced in the code above. Use the status column (instead of binary_status) which has censored, vte, death
fg_data <- finegray(Surv(time_to_status_months, status) ~ .,
                    data = fit_data, 
                    etype = 1,
                    id = participant_id)
# run the cox model. combined category is the PRS/PCDH15 combined category
model <- coxph(Surv(fgstart, fgstop, fgstatus) ~combined_category,
                      data = fg_data)

# create the new_df for surv-fit
new_df <- data.frame(combined_category = as.factor(c('high PRS; PCDH15 mutated',
                                                     'high PRS; PCDH15 wiltype',
                                                     'low PRS; PCDH15 mutated',
                                                     'low PRS; PCDH15 wildtype'
)))
fit <- survfit(model, newdata = new_df)
# invert the fit to get cuminc instead of survival
fit$surv=1-fit$surv
# for more flexible altering
surv_data <- surv_summary(fit)
# map the correct term to each strata
levels(surv_data$strata) <- c('high PRS; PCDH15 mutated',
                              'high PRS; PCDH15 wiltype',
                              'low PRS; PCDH15 mutated',
                              'low PRS; PCDH15 wildtype')

write.table(surv_data, './results/cuminc_plot_PRS_PCDH15_cmprsk.tsv',
            row.names = F, quote = F, sep = '\t')


# TMB interactions #####
# load the file with the TMBs and mutational signatures 
tmb_ac <- fread('./data/tmb_and_mutsig_ac.csv') %>% 
  # select only the columns relating to TMB
  select(tumour_sample_platekey, somatic_coding_variants_per_mb, starts_with('tmb'), starts_with('signature')) %>%
  # create a column tmb_binary which is high if somatic_coding_mutations_per_mb >= 10
  mutate(tmb_binary = factor(case_when(somatic_coding_variants_per_mb >= 20 ~ 'high_thresh20', .default = 'low'), levels = c('low', 'high_thresh20')))


temp_data <- left_join(tdata, tmb_ac, by = 'tumour_sample_platekey')

# model with TMB and PRS
results <- coxph(Surv(time = time_to_status, event = binary_status) ~ tmb_binary + klarin_PRS_scaled + tmb_binary:klarin_PRS_scaled +
                          age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                        data = temp_data) %>% 
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
    filter(str_detect(term, 'tmb_binary')|str_detect(term, 'klarin')) %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, HR, l95ci, u95ci, pval) %>%
    mutate(model_covars = ' tmb_binary + klarin_PRS_scaled + tmb_binary:klarin_PRS_scaled + age_study_entry + sex + pc1-4')

write.table(results, './results/cox_PRS_TMB_interactions.tsv', row.names=F, quote = F, sep = '\t')

# signature interactions ####
tmb_ac <- fread('./data/tmb_and_mutsig_ac.csv') %>% 
  # select only the columns relating to signatures
  select(tumour_sample_platekey, starts_with('signature')) %>%
  # NA means signature contributed <5% to TMB: then code as zero, otherwise code as 1
  mutate(across(starts_with('signature'), ~ifelse(is.na(.), 0,1)))

temp_data <- tdata %>% 
  # merge with the TMB data
  left_join(., tmb_ac, by = 'tumour_sample_platekey') %>%
  # add +0.5 to time_to_status only where time_to_status =0
  # this step is required becuse tmerge does not accept stop times of 0
  # see vignette by Therneau et al
  mutate(time_to_status = case_when(time_to_status == 0 ~ 0.5,
                                    .default = time_to_status ))

results <- list()
for (i in c('signature_8', 'signature_6', 'signature_19', 'signature_26')){
results[[i]] <- coxph(as.formula(paste0('Surv(time = time_to_status, event = binary_status) ~ ', 
                                       i, ' + klarin_PRS_scaled +', i, ':klarin_PRS_scaled +
                   age_study_entry + sex + pc1 + pc2 + pc3 + pc4')), 
                 data = temp_data) %>% 
  # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
  tidy(exp=T, conf.int = T) %>% 
  # round estimates to 4 decimals
  mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
  # the output provides the HR for VTE associated with each of the exposure variables in the formula,
  # after adjusting for all the other variables
  # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
  filter(str_detect(term, i)|str_detect(term, 'klarin')) %>% 
  # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
  rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
  # select relevant columns for output
  select(term, HR, l95ci, u95ci, pval) %>%
  mutate(model_covars = paste0(i,' + klarin_PRS_scaled + ',i, ':klarin_PRS_scaled + age_study_entry + sex + pc1-4'))
}

results_df <- do.call(rbind, results)
write.table(results_df, './results/cox_PRS_signature_interactions.tsv', row.names=F, quote = F, sep = '\t')

# END
