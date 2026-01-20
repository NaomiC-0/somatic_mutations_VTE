# Regression models
# N Cornish
# script for running competing risk regression models with cancer data 

# Config ####
# load config script
rm(list=ls())
setwd(Sys.getenv("wd"))
source('./scripts/config421.R')
library(cmprsk)
library(broom)
library(survminer)

# Load formatted analysis cohort ####
ac <- readRDS('./data/analysis_cohort_formatted.rds')

table(ac$status)
table(ac$status_description)

# Li et al risk classification is as follows
# +3 Pancreatic, gastric, esophageal, cholangiocarcinoma, and gallbladder 
# +2 Lung, ovarian, uterine, bladder, kidney, testicular, aggressive NHL, myeloma, brain, and soft tissue sarcoma 
# +1 Colorectal/intestinal cancer 
# 0 : all other


#Create sdata frame:
# create sdata where sact is coded as a time varying covariate
#tmerge 
# use tmerge because sact is a time dependent covar
# first create the basic dataset for tmerge
tdata <- ac %>% 
  mutate(li_risk_group = case_when(
    # +3 VERY HIGH RISK
    site_coded_3char %in% c(
      "C25",  # Pancreas
      "C15",  # Oesophagus
      "C16",  # Stomach
      "C23",  # Gallbladder
      "C24"   # Extrahepatic bile ducts / cholangiocarcinoma
    ) ~ "Very high risk", 
    
    # +2 HIGH RISK
    site_coded_3char %in% c(
      "C34",  # Lung
      "C56",  # Ovary
      "C54",  # Uterus (endometrium)
      "C55",  # Uterus NOS (cervix excluded)
      "C67",  # Bladder
      "C64",  # Kidney
      "C62",  # Testis
      "C83",  # Aggressive non-Hodgkin lymphoma
      "C85",  # NHL NOS
      "C90",  # Multiple myeloma
      "C71",  # Brain
      "C49",  # Soft tissue sarcoma
      "C40",  # Bone sarcoma (long bones)
      "C41",  # Bone sarcoma (other bones)
      "C42"   # Haematological malignancy
    ) ~ "High risk",
    
    # +1 INTERMEDIATE RISK
    site_coded_3char %in% c(
      "C18",  # Colon
      "C19",  # Rectosigmoid
      "C20",  # Rectum
      "C17"  # Small intestine
    ) ~ "Intermediate risk",
    # 0 LOW RISK
    .default = "Low risk"
  )) %>%
  mutate(li_risk_group = factor(li_risk_group, levels = c("Low risk",
                                "Intermediate risk",
                                "High risk",
                                "Very high risk" )),
# add +0.5 to time_to_status only where time_to_status =0
  # this step is required becuse tmerge does not accept stop times of 0
  # see vignette by Therneau et al
  time_to_status = case_when(time_to_status == 0 ~ 0.5,
                                    .default = time_to_status ))
# now perform tmerge to split f/u times for participants who received sact after study entry
sdata <- tmerge(data1 = tdata, data2 = tdata,
                vte = event(time_to_status, status), 
                current_sact = tdc(time_to_current_sact), 
                id = participant_id,
                options= list(idname="participant_id")) %>%
  # select only relevant columns
  select(participant_id, tumour_sample_platekey, age_study_entry,
         sex, pc1, pc2, pc3, pc4,
         disease_type, stage_grouped,
         sact_over6weeks_before_studyentry, current_sact,
         tstart, tstop, vte)

# Identify only genes which were possibly assoc with VTE in the cox models
gene_list <- fread('./results/cox_genes_minadj.tsv') %>%
  # select genes with pval <= 0.05 in primary analysis
  filter(fdr_p<=0.1) %>% distinct(gene) %>% unlist

# Load somatic variants ####
all_variants <-fread('./data/combined_variants.csv') %>% 
  # participants will only be counted once for each gene even if they carry multiple different variant
  distinct

# Load germline PRS
prs_scores <- fread('./data/participant_PRS_scores_klarin.txt.sscore') %>%
  # filter only for the IID column and SCORE1_AVG column. Note the average means the total PRS was then divided by the number of SNPs - this denominator was 299 if a participant had a valid genotype for every single SNP but due to quality control most participants had genotypes missing for a number of the SNPs.
  dplyr::select(IID, SCORE1_AVG) %>%
  # change column names so they are more intuitive
  dplyr::rename('germline_sample_platekey' = IID,
                'klarin_PRS' = SCORE1_AVG) %>%
  filter(germline_sample_platekey %in% ac$germline_sample_platekey) %>%
  # scale the PRS to have a mean of zero and SD of 1
  mutate(klarin_PRS_scaled = scale(klarin_PRS) %>% as.vector) %>%
  # when you scale need to convert column class back to vector as otherwise your column names appear as colname[,1]
  # create a PRS quartile column
  mutate(prs_quartile = ntile(klarin_PRS_scaled, 4)) %>%
  mutate(prs_category = case_when(prs_quartile == 1 ~ 'lower_3_quartiles',
                                  
                                  prs_quartile == 2 ~ 'lower_3_quartiles',
                                  
                                  prs_quartile == 3 ~ 'lower_3_quartiles',
                                  
                                  prs_quartile == 4 ~ 'top_quartile',
                                  
                                  .default = NA
                                  
  )) %>%
  select(germline_sample_platekey, prs_category)

# Genes FG minadj ####
results <- list()
for (i in gene_list) {
  # filter the all_variants table for the gene in question
  variants_by_gene <- all_variants %>% filter(gene == i)
  # create a new column in the analysis cohort dataframe 
  temp_data <- sdata %>% 
    # which identifies whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  # fine gray weights
  fg_data <- finegray(Surv(tstart, tstop, vte) ~ gene + age_study_entry + sex + pc1 + pc2 + pc3 + pc4,
                      data = temp_data, 
                      etype = 1,
                      id = participant_id)
    results[[i]] <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 
                        gene + age_study_entry + sex + pc1 + pc2 + pc3 + pc4,
                      data = fg_data) %>% 
    tidy(exp = T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
    filter(term == 'gene') %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    dplyr::rename('gene' = term, 'SHR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(gene, SHR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(gene = i,
           model_covars = paste0('fine_gray(time,status) ~ ', i, ' + age_study_entry + sex + pc1-4'))}

results_df <- do.call(rbind, results) 
write.table(results_df, './results/FG_genes_minadj.tsv', row.names=F, quote = F, sep = '\t')

# Gene FG fulladj ####
results <- list()
for (i in gene_list) {
  # filter the all_variants table for the gene in question
  variants_by_gene <- all_variants %>% filter(gene == i)
  # create a new column in the analysis cohort dataframe 
  temp_data <- sdata %>% 
    # which identifies whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  fg_data <- finegray(Surv(tstart, tstop, vte) ~ .,
                      data = temp_data, 
                      etype = 1,
                      id = participant_id)
  results[[i]] <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 
                          gene + age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                          disease_type + stage_grouped + 
                          sact_over6weeks_before_studyentry + current_sact, 
                        data = fg_data) %>% 
    tidy(exp = T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
    filter(term == 'gene') %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    dplyr::rename('gene' = term, 'SHR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(gene, SHR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(gene = i,
           model_covars = paste0('fine_gray(time,status) ~ ', i, ' + age_study_entry + sex + pc1-4 + disease_type + stage_grouped + sact_over6weeks_before_studyentry + current_sact + '))}

results_df <- do.call(rbind, results) 
write.table(results_df, './results/FG_genes_fulladj.tsv', row.names=F, quote = F, sep = '\t')

# Cuminc whole cohort #####
fg_data <- finegray(Surv(time_to_status_months, status) ~ .,
                    data = tdata, 
                    etype = 1,
                    id = participant_id)
model <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1, data = fg_data)

fit <- survfit(model)
fit2 <- fit
# invert the fit to get cuminc instead of survival
fit2$surv=1-fit$surv
fit2$upper=1-fit$lower
fit2$lower=1-fit$upper
# for more flexible altering
surv_data <- surv_summary(fit2)

write.table(surv_data, './results/cuminc_plot_whole_cohort.tsv', row.names = F,
            quote = F, sep = '\t')

# Genes Cuminc #####
# for the cuminc curves model by month
create_surv_dat <- function(gene_name){
  variants_by_gene <- all_variants %>% filter(gene == gene_name)
  # create a new column in the analysis cohort dataframe 
  temp_data <- ac %>%
    # which identifies whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  fg_data <- finegray(Surv(time_to_status_months, status) ~ .,
                      data = temp_data, 
                      etype = 1,
                      id = participant_id)
  model <- coxph(Surv(fgstart, fgstop, fgstatus) ~ gene, data = fg_data)
  new_df <- data.frame(gene = as.factor(c(0,1)))
  fit <- survfit(model, newdata = new_df)
  # invert the fit to get cuminc instead of survival
  fit2$surv=1-fit$surv
  # for more flexible altering
  surv_data <- surv_summary(fit2)
  # map the correct term to each strata
  levels(surv_data$strata) <- c('wildtype','mutated')
  return(surv_data)
  }

for (i in gene_list) {
 df <- create_surv_dat(i)
 write.table(df, paste0('cuminc_plot_', i, '.tsv'),
           row.names = F, quote = F, sep = '\t')
}

# Li-Cancer risk group cuminc plot ######
# apply fine gray weights for competing risks
fg_dat <- finegray(
  Surv(time_to_status_months, status) ~ li_risk_group,
  data = tdata,
  etype = 1
)
# fit cox model
model <- coxph(
  Surv(fgstart, fgstop, fgstatus) ~ li_risk_group,
  data = fg_dat
)

## cindex list ####
cindex <- list()
cindex[["tumour_type"]] <- data.frame(
  cindex = (concordance(model))$concordance,
  model = 'Li cancer risk category only')

# 
new_df <- data.frame(
  li_risk_group = as.factor(c("Low risk",
                              "Intermediate risk",
                              "High risk",
                              "Very high risk" 
)))
# get cumincs
fit <- survfit(model, newdata = new_df)
# invert the fit to get cuminc instead of survival
fit2 <- fit
fit2$surv=1-fit$surv
fit2$lower=1-fit$upper
fit2$upper=1-fit$lower
# for more flexible altering
surv_data <- surv_summary(fit2)
# map the correct term to each strata
levels(surv_data$strata) <- c("Low risk",
                              "Intermediate risk",
                              "High risk",
                              "Very high risk" 
)

write.table(surv_data, paste0('./results/cuminc_plot_Li_cancer_grps_cmprsk.tsv'),
            row.names = F, quote = F, sep = '\t')

# number in each cancer risk group
table(tdata$li_risk_group) %>% as.data.frame() %>%
  write.table(., './results/freqs_li_cancer_risk_grp.tsv', row.names = F, sep = '\t',
              quote =F)
rm(fg_dat, fit, model, new_df, surv_data)

# Genes: Li-cancer category interactions #####
create_surv_dat <- function(gene_name){
variants_by_gene <- all_variants %>% filter(gene == gene_name)
# create a new column in the analysis cohort dataframe 
temp_data <- tdata %>%
  # which identifies whether the participants sample appears in the variants_by_gene table
  mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0),
         ) %>%
  # create a new category which combines Kscore status and somatic mutation status
  mutate(combined_category = case_when(
    # category 1
    li_risk_group == 'Very high risk' & gene == '1' ~ 
      'very_high_risk_cancer_plus_somatic_mutation',
    # category 2 
    li_risk_group == 'Very high risk' & gene == '0' ~ 
      'very_high_risk_cancer_no_mutation',
    # category 3
    li_risk_group == 'High risk' & gene == '1' ~ 
      'high_risk_cancer_plus_somatic_mutation',
    # category 4
    li_risk_group == 'High risk' & gene == '0' ~ 
      'high_risk_cancer_no_mutation',
    # category 5
    li_risk_group == 'Intermediate risk' & gene == '1' ~ 
      'intermediate_risk_cancer_plus_somatic_mutation',
    # category 6
    li_risk_group == 'Intermediate risk' & gene == '0' ~ 
      'intermediate_risk_cancer_no_mutation',
    # category 7
    li_risk_group == 'Low risk' & gene == '1' ~ 
      'low_risk_cancer_plus_somatic_mutation',
    # category 8
    li_risk_group == 'Low risk' & gene == '0' ~ 
      'low_risk_cancer_no_mutation'),
 # create factor levels
    combined_category = factor(combined_category, levels = c(
      "very_high_risk_cancer_plus_somatic_mutation",
      "very_high_risk_cancer_no_mutation",
      "high_risk_cancer_plus_somatic_mutation",
      "high_risk_cancer_no_mutation",
      "intermediate_risk_cancer_plus_somatic_mutation",
      "intermediate_risk_cancer_no_mutation",
      "low_risk_cancer_plus_somatic_mutation",
      "low_risk_cancer_no_mutation")
    ))

# create table of the numbers of participants in each group and write to file
table(temp_data$combined_category) %>% as.data.frame() %>%
  write.table(., paste0('./results/freqs_li_cancer_risk_grp_plus', gene_name, '.tsv'), row.names = F, sep = '\t',
              quote =F)
fg_data <- finegray(Surv(time_to_status_months, status) ~ .,
                    data = temp_data, 
                    etype = 1,
                    id = participant_id)
# run the cox model
model <- coxph(Surv(fgstart, fgstop, fgstatus) ~ combined_category,
               data = fg_data)

# create the new_df for surv-fit
new_df <- data.frame(combined_category = as.factor(c(
                                                     "very_high_risk_cancer_plus_somatic_mutation",
                                                     "very_high_risk_cancer_no_mutation",
                                                     "high_risk_cancer_plus_somatic_mutation",
                                                     "high_risk_cancer_no_mutation",
                                                     "intermediate_risk_cancer_plus_somatic_mutation",
                                                     "intermediate_risk_cancer_no_mutation",
                                                     "low_risk_cancer_plus_somatic_mutation",
                                                     "low_risk_cancer_no_mutation"
                                                     )))
fit <- survfit(model, newdata = new_df)

# invert the fit to get cuminc instead of survival
fit2 <- fit
fit2$surv=1-fit$surv
fit2$upper=1-fit$lower
fit2$lower=1-fit$upper
# for more flexible altering
surv_data <- surv_summary(fit2)
# map the correct term to each strata
levels(surv_data$strata) <- c("very_high_risk_cancer_plus_somatic_mutation",
                              "very_high_risk_cancer_no_mutation",
                              "high_risk_cancer_plus_somatic_mutation",
                              "high_risk_cancer_no_mutation",
                              "intermediate_risk_cancer_plus_somatic_mutation",
                              "intermediate_risk_cancer_no_mutation",
                              "low_risk_cancer_plus_somatic_mutation",
                              "low_risk_cancer_no_mutation"
)

# save surv_fit table
write.table(surv_data, paste0('./results/cuminc_plot_Li_cancer_category_', gene_name, '_cmprsk.tsv'),
            row.names = F, quote = F, sep = '\t')

# return cindex for model
cindex[[gene_name]] <<- data.frame(
  cindex = (concordance(model))$concordance,
  model = paste(gene_name, 'mutation status + Li cancer risk category')
)
}

for (gene_name in c('CDKN2A', 'KRAS', 'PCDH15', 'TP53')) {
  create_surv_dat(gene_name)
  }
cindex_df <- do.call(rbind, cindex)
write.table(cindex_df, './results/cindex_li_cancer_category_plus_somatic_mutation.tsv',
            quote = F, sep = '\t')

# Stage cuminc plot ######
#  exclude patients with unknown stage
# apply fine gray weights for competing risks
tdata_stage_model <- tdata %>% filter(stage_grouped != 'unknown')

fg_dat <- finegray(
  Surv(time_to_status_months, status) ~ stage_grouped,
  data = tdata_stage_model,
  etype = 1
)
# fit cox model
model <- coxph(
  Surv(fgstart, fgstop, fgstatus) ~ stage_grouped,
  data = fg_dat
)

## cindex list ####
cindex <- list()
cindex[["stage_grouped"]] <- data.frame(
  cindex = (concordance(model))$concordance,
  model = 'stage only')

# 
new_df <- data.frame(
  stage_grouped = as.factor(c("early","late" 
  )))
# get cumincs
fit <- survfit(model, newdata = new_df)
# invert the fit to get cuminc instead of survival
fit2$surv=1-fit$surv
fit2$lower=1-fit$upper
fit2$upper=1-fit$lower
# for more flexible altering
surv_data <- surv_summary(fit2)
# map the correct term to each strata
levels(surv_data$strata) <- c("early","late")

write.table(surv_data, paste0('./results/cuminc_plot_stage_only_cmprsk.tsv'),
            row.names = F, quote = F, sep = '\t')

# number in each cancer risk group
table(tdata_stage_model$stage_grouped) %>% as.data.frame() %>%
  write.table(., './results/freqs_stage_grp.tsv', row.names = F, sep = '\t',
              quote =F)
rm(fg_dat, fit, model, new_df, surv_data)

# Genes:stage #####
# now examine interactions
create_surv_dat_stage <- function(gene_name){
  variants_by_gene <- all_variants %>% filter(gene == gene_name)
  # create a new column in the analysis cohort dataframe 
  temp_data <- tdata_stage_model %>%
    # which identifies whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0),
    ) %>%
    # create a new category which combines Kscore status and somatic mutation status
    mutate(combined_category = case_when(stage_grouped == 'early' &
                                           gene == '1' ~ 'early_stage_cancer_plus_somatic_mutation',
                                         stage_grouped == 'early' &
                                           gene == '0' ~ 'early_stage_cancer',
                                         stage_grouped == 'late' &
                                           gene == '1' ~ 'late_stage_cancer_plus_somatic_mutation',
                                         stage_grouped == 'late' &
                                           gene == '0' ~ 'late_stage_cancer',
                                 #        stage_grouped == 'unknown' &
                                  #         gene == '1' ~ 'unknown_stage_cancer_plus_somatic_mutation',
                                   #      stage_grouped == 'unknown' &
                                    #       gene == '0' ~ 'unknown_stage_cancer', 
                                 .default = NA),
           combined_category = factor(combined_category, levels = c('late_stage_cancer_plus_somatic_mutation',
                                                                    'late_stage_cancer',
                                                         #           'unknown_stage_cancer_plus_somatic_mutation',
                                                          #          'unknown_stage_cancer',
                                                                    'early_stage_cancer_plus_somatic_mutation',
                                                                    'early_stage_cancer'
           )))
  # frequency table
  table(temp_data$combined_category) %>% as.data.frame() %>%
    write.table(., paste0('./results/freqs_stage_plus', gene_name, '.tsv'), row.names = F, sep = '\t',
                quote =F)
  
  fg_data <- finegray(Surv(time_to_status_months, status) ~ combined_category,
                      data = temp_data, 
                      etype = 1,
                      id = participant_id)
  # run the cox model
  model <- coxph(Surv(fgstart, fgstop, fgstatus) ~ combined_category,
                 data = fg_data)
  # return cindex for model
  cindex[[gene_name]] <<- data.frame(
    cindex = (concordance(model))$concordance,
    model = paste(gene_name, 'mutation status + cancer stage')
  )
  # create the new_df for surv-fit
  new_df <- data.frame(combined_category = as.factor(c('late_stage_cancer_plus_somatic_mutation',
                                                       'late_stage_cancer',
                                              #         'unknown_stage_cancer_plus_somatic_mutation',
                                              #         'unknown_stage_cancer',
                                                       'early_stage_cancer_plus_somatic_mutation',
                                                       'early_stage_cancer'
  )))
  fit <- survfit(model, newdata = new_df)
  # invert the fit to get cuminc instead of survival
  fit2 <- fit
  fit2$surv=1-fit$surv
  fit2$lower=1-fit$upper
  fit2$upper=1-fit$lower
  # for more flexible altering
  surv_data <- surv_summary(fit2)
  # map the correct term to each strata
  levels(surv_data$strata) <- c('late_stage_cancer_plus_somatic_mutation',
                                'late_stage_cancer',
                          #      'unknown_stage_cancer_plus_somatic_mutation',
                          #      'unknown_stage_cancer',
                                'early_stage_cancer_plus_somatic_mutation',
                                'early_stage_cancer'
  )
  
  write.table(surv_data, paste0('./results/cuminc_plot_stage_', gene_name, '_cmprsk.tsv'),
              row.names = F, quote = F, sep = '\t')
}
 
for (i in c('TP53', 'CDKN2A', 'KRAS', 'PCDH15')) {
  create_surv_dat_stage(i)
} 

cindex_df <- do.call(rbind, cindex)
write.table(cindex_df, './results/cindex_stage_plus_somatic_mutation.tsv',
            quote = F, sep = '\t')


# TMB ####
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
    levels = c('tmb_below_5','tmb_between_5_10','tmb_between_10_20', 'tmb_above_20')),
  )
# tmb list
tmb_list <- c('tmb_binary_thresh20', 'tmb_categories')

temp_data <- left_join(sdata, tmb_ac, by = join_by("tumour_sample_platekey"))
## minadj ####
results <- list()
for (i in tmb_list) {
  fg_data <- finegray(Surv(tstart, tstop, vte) ~ .,
                      data = temp_data, 
                      etype = 1,
                      id = participant_id)
  results[[i]] <- coxph(as.formula(paste('Surv(fgstart, fgstop, fgstatus) ~', 
                                         i,  '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4')),
                        data = fg_data) %>% 
    tidy(exp = T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
    filter(str_detect(term, i)) %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    dplyr::rename('SHR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, SHR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(model_covars = paste0('fine_gray(time,status) ~ ', i, ' + age_study_entry + sex + pc1-4'))}

results_df <- do.call(rbind, results) 
write.table(results_df, './results/FG_TMB_minadj.tsv', row.names=F, quote = F, sep = '\t')

## fulladj ####
results <- list()
for (i in tmb_list) {
  fg_data <- finegray(Surv(tstart, tstop, vte) ~ .,
                      data = temp_data, 
                      etype = 1,
                      id = participant_id)
  results[[i]] <- coxph(as.formula(paste('Surv(fgstart, fgstop, fgstatus) ~', 
                                         i,  '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + disease_type + stage_grouped + 
                          sact_over6weeks_before_studyentry + current_sact')),
                        data = fg_data) %>% 
    tidy(exp = T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
    filter(str_detect(term, i)) %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    dplyr::rename('SHR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, SHR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(model_covars = paste0('fine_gray(time,status) ~ ', i, ' + age_study_entry + sex + pc1-4'))}

results_df <- do.call(rbind, results) 
write.table(results_df, './results/FG_TMB_fulladj.tsv', row.names=F, quote = F, sep = '\t')

# signatures ####
# load the file with the TMBs and mutational signatures 
tmb_ac_continuous <- fread('./data/tmb_and_mutsig_ac.csv') %>% 
  # select only the columns relating to signatures
  select(tumour_sample_platekey, starts_with('signature'))
# convert signatures (continuous numerical with high proprtion NA/zero vals to a binary val)
# if NA then code as zero, otherwise code as 1
tmb_ac <- tmb_ac_continuous %>% mutate(across(starts_with('signature'), ~ifelse(is.na(.), 0,1)))

# sig list
sig_list <- c('signature_6', 'signature_8', 'signature_19', 'signature_26')

temp_data <- left_join(sdata, tmb_ac, by = join_by("tumour_sample_platekey"))
## minadj ####
results <- list()
for (i in sig_list) {
   fg_data <- finegray(Surv(tstart, tstop, vte) ~ .,
                      data = temp_data, 
                      etype = 1,
                      id = participant_id)
  results[[i]] <- coxph(as.formula(paste('Surv(fgstart, fgstop, fgstatus) ~', 
                        i,  '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4')),
                        data = fg_data) %>% 
    tidy(exp = T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
    filter(term == i) %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    dplyr::rename('SHR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, SHR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(model_covars = paste0('fine_gray(time,status) ~ ', i, ' + age_study_entry + sex + pc1-4'))}

results_df <- do.call(rbind, results) 
write.table(results_df, './results/FG_signatures_minadj.tsv', row.names=F, quote = F, sep = '\t')

## fulladj ####
results <- list()
for (i in sig_list) {
  fg_data <- finegray(Surv(tstart, tstop, vte) ~ .,
                      data = temp_data, 
                      etype = 1,
                      id = participant_id)
  results[[i]] <- coxph(as.formula(paste('Surv(fgstart, fgstop, fgstatus) ~', 
                                         i,  '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4 + disease_type + stage_grouped + 
                          sact_over6weeks_before_studyentry + current_sact')),
                        data = fg_data) %>% 
    tidy(exp = T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
    filter(term == i) %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    dplyr::rename('SHR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, SHR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(model_covars = paste0('fine_gray(time,status) ~ ', i, 
                                 ' + age_study_entry + sex + pc1-4 + disease_type + stage_grouped + sact_over6weeks_before_studyentry + current_sact'))}

results_df <- do.call(rbind, results) 
write.table(results_df, './results/FG_signatures_fulladj.tsv', row.names=F, quote = F, sep = '\t')



