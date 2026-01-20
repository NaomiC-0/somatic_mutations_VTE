# Cox regression models for gene level analyses
# N Cornish

# Config ####
# load config script
setwd(Sys.getenv("wd"))
source('./scripts/config421.R')

# Load formatted analysis cohort ####
ac <- readRDS('./data/analysis_cohort_formatted.rds')

# Load the gene_list_for_analysis produced using the script gene_list.R
gene_list <- fread('./data/gene_list_for_analysis.txt') %>% 
  distinct(gene) %>% unlist
# Load genetic variants ####
all_variants <-fread('./data/combined_variants.csv') %>% 
  # participants will only be counted once for each gene even if they carry multiple different variant
  distinct

# Cox models ####
## minimally adjusted model #####
# adjusted for age, sex and genetic PCs 1- 4

# create list to store results
results <- list()
#loop through the gene list to run the following function on each gene individually
for (i in gene_list) {
  # filter the all_variants table for the gene in question
  variants_by_gene <- all_variants %>% filter(gene == i)
  # create a new column in the analysis cohort dataframe 
  temp_data <- ac %>% 
    # which identifies whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  # now run the cox model
  results[[i]] <- coxph(Surv(time = time_to_status, event = binary_status) ~ gene + age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                        data = temp_data) %>% 
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
    filter(term == 'gene') %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    rename('gene' = term, 'HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(gene, HR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(gene = i)
}

# merge into a single df and look at the fdr adjusted p values
results_df <- do.call(rbind, results) %>%
  # add adjusted p values
  mutate(fdr_p = p.adjust(pval, method = 'fdr')) %>% 
  arrange(fdr_p) %>%
  # add a column which states the model run
  mutate(model_covars = 'gene + age_study_entry + sex + pc1-4')

write.table(results_df, './results/cox_genes_minadj.tsv', row.names=F, quote = F, sep = '\t')

## Fully adjusted model ####

#Create sdata frame:
# create sdata where sact is coded as a time varying covariate
#tmerge 
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
  select(participant_id, tumour_sample_platekey, age_study_entry,
         sex, pc1, pc2, pc3, pc4,
         disease_type, stage_grouped,
         sact_over6weeks_before_studyentry, current_sact,
         tstart, tstop, vte)

results <- list()
#loop through the gene list to run the following function on each gene individually
for (i in gene_list) {
  # filter the all_variants table for the gene in question
  variants_by_gene <- all_variants %>% filter(gene == i)
  # create a new column in the analysis cohort dataframe NB using sdata here, 
  # where participants who receive sact after study entry have split time intervals
  temp_data <- sdata %>% 
    #  identifies whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  # now run the cox model
  results[[i]] <- coxph(Surv(tstart, tstop, vte) ~ gene + 
                          age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                          disease_type + stage_grouped + 
                          sact_over6weeks_before_studyentry + current_sact, 
                        data = temp_data) %>%
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4 and disease_type)
    filter(term == 'gene') %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    rename('gene' = term, 'HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(gene, HR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(gene = i)
}

# merge into a single df and look at the fdr adjusted p values
results_df <- do.call(rbind, results) %>%
  arrange(pval) %>%
  # add a column which states the model run
  mutate(model_covars = 'gene + age_study_entry + sex + pc1-4 + disease_type + stage_grouped + sact_over6weeks_before_studyentry + current_sact')


write.table(results_df, './results/cox_genes_fulladj.tsv', row.names=F, quote = F, sep = '\t')

## multivariable gene model ####
#for the genes with fdr_p < 0.1 check whether the effects are independent of each other in an additive model

## fdr_sig_genes #####
fdr_sig_genes <- fread('./results/cox_genes_minadj.tsv') %>% filter(fdr_p <= 0.1) %>% pull(gene)
print(fdr_sig_genes)

# create a new function which create a column with the genotypes for all the genes in fdr_sig_genes
genotypes <- function(fdr_sig_genes){
  df <- ac
  #loop through the gene list to run the following function on each gene individually
  for (i in fdr_sig_genes) {
    # filter the all_variants table for the gene in question
    variants_by_gene <- all_variants %>% filter(gene == i)
    # create a new column in the analysis cohort dataframe 
    df[[i]] <- ifelse(df$tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0)
  }
  return(df)}
# create temp_data with the genotyeps for the fdr_sig genes
temp_data <- genotypes(fdr_sig_genes)  
# check
table(temp_data$TP53)
table(temp_data$KRAS)


# now run the cox model
results <- coxph(Surv(time = time_to_status, event = binary_status) ~ TP53 + PCDH15 + KRAS + CDKN2A,
                 data = temp_data) %>% 
  # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
  tidy(exp=T, conf.int = T) %>% 
  # round estimates to 4 decimals
  mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
  # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
  rename('gene' = term, 'HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
  # select relevant columns for output
  select(gene, HR, l95ci, u95ci, pval) %>%
  # add a column which states the model run
  mutate(model_covars = 'TP53 + PCDH15 + KRAS + CDKN2A')

write.table(results, './results/cox_multigene_regression.tsv', row.names=F, quote = F, sep = '\t')


# Tumour stratified analysis ####
# run this tumour stratified analysis only on fdr_sig_genes

# first create vector of disease_types
tumour_types <- ac %>% 
  # remove 'OTHER' because insufficient VTE cases to run stratified analysis 
  filter(disease_type != 'OTHER') %>%
  # create vector of remaining disease types
  pull(disease_type) %>% unique

results <- list()
# for each disease type run the analysis
for (t in tumour_types) {
  # filter the analysis cohort for tumour type t
  sdata_tumour <- sdata %>% filter(disease_type == t)
  # create a new list to store results for that tumour
  results[[t]] <- list()
  # now loop through the gene list to run the following function on each gene individually
  for (i in fdr_sig_genes) {
    # filter the all_variants table for the gene in question
    variants_by_gene <- all_variants %>% filter(gene == i)
    # create a new column in the temp_data frame
    temp_data <- sdata_tumour %>% 
      #  identifies whether the participants sample appears in the variants_by_gene table
      mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
    # now run the cox model adjusting for everything except tumour type
    results[[t]][[i]] <- coxph(Surv(tstart, tstop, vte) ~ gene + 
                                 age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                                 stage_grouped + 
                                 sact_over6weeks_before_studyentry + current_sact, 
                               data = temp_data) %>%
      # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
      tidy(exp=T, conf.int = T) %>% 
      # round estimates to 4 decimals
      mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
      # the output provides the HR for VTE associated with each of the exposure variables in the formula,
      # after adjusting for all the other variables
      # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4 and disease_type)
      filter(term == 'gene') %>% 
      # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
      rename('gene' = term, 'HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
      # select relevant columns for output
      select(gene, HR, l95ci, u95ci, pval) %>%
      # replace 'gene' with the actual name of the gene queried
      mutate(gene = i, cancer = t)
  }}

# note had issues with model convergence for some tumour-gene combinations due to data sparsity
results_df <- lapply(results, function(x) do.call(rbind, x)) %>%
  do.call(rbind,.) %>%
  arrange(cancer, pval) %>%
  # add a column which states the model run
  mutate(model_covars = 'gene + age_study_entry + sex + pc1-4 + stage_grouped + sact_over6weeks_before_studyentry + current_sact') %>%
  # a few of the cancers have uci95 = Inf as the model did not converge due to sparse data
  # for these cancers the analysis is invalid therefore substitute with NA vals
  mutate(HR = case_when(u95ci == 'Inf' ~ NA, .default = HR),
         l95ci = case_when(u95ci == 'Inf' ~ NA, .default = l95ci),
         pval = case_when(u95ci == 'Inf' ~ NA, .default = pval),
         u95ci = case_when(u95ci == 'Inf' ~ NA, .default = u95ci)
  )


write.table(results_df, './results/cox_genes_tumour_stratified.tsv', row.names=F, quote = F, sep = '\t')

# Tumour interactions #### 
##  ANOVA #####
results_anova <- list()
#loop through the gene list to run the following function on each gene individually
for (i in c('TP53', 'KRAS', 'CDKN2A', 'PCDH15', 'IDH1')) {
  # filter the all_variants table for the gene in question
  variants_by_gene <- all_variants %>% filter(gene == i)
  # create a new column in the analysis cohort dataframe NB using sdata here, 
  # where participants who receive sact after study entry have split time intervals
  temp_data <- sdata %>% 
    #  identifies whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  # additive model
  additive_model <- coxph(Surv(tstart, tstop, vte) ~ gene + disease_type +
                            age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                            stage_grouped + 
                            sact_over6weeks_before_studyentry + current_sact, 
                          data = temp_data)
  # interaction cox model
  interaction_model <- coxph(Surv(tstart, tstop, vte) ~ gene + 
                               disease_type + disease_type:gene +
                          age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                          stage_grouped + 
                          sact_over6weeks_before_studyentry + current_sact, 
                        data = temp_data) 
  results_anova[[i]] <- anova(additive_model, interaction_model)  
}
  # extract anova pvalue resutls to table and save
results_anova_df <- tibble(
  gene = names(results_anova),
  p_value = sapply(results_anova, function(x) x$`P(>|Chi|)`[2])
) %>%
  mutate(test = paste('anova(additive model with disease_type +', 
                      gene, 'vs interaction model incorporating disease_type :', gene))

write.table(results_anova_df, './results/anova_tumour_gene_interactions.tsv',
            row.names = F, quote = F, sep = '\t')
rm(results_anova_df, results_anova, interaction_model, additive_model)

## INDIVIDUAL ####
# NEXT do the same thing but extract interaction terms for EACH tumour type, 
# BUT COMPARE with the REFERENCE of 'all OTHER' tumour types 
#because otherwise the individual interaction term pvals fluctuate a lot depending on which tumour is used as ref
results <- list()
for (i in c('TP53', 'KRAS', 'CDKN2A', 'PCDH15')) {
  # filter the all_variants table for the gene in question
  variants_by_gene <- all_variants %>% filter(gene == i)
  # create a new column in the analysis cohort dataframe NB using sdata here, 
  # where participants who receive sact after study entry have split time intervals
  temp_data <- sdata %>% 
    #  identifies whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  for (t in unique(ac$disease_type)) {
  # make all other disease groups the reference
    temp_data_tumour_binary <- temp_data %>%
      mutate(cancer_type_binary = case_when(disease_type == t ~ 1, .default = 0))
    # interaction cox model
  results[[i]][[t]] <- coxph(Surv(tstart, tstop, vte) ~ gene + 
                               cancer_type_binary +
                               cancer_type_binary:gene +
                               age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                               + stage_grouped + 
                               sact_over6weeks_before_studyentry + current_sact, 
                             data = temp_data_tumour_binary) %>%
# use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4 and disease_type)
    filter(str_detect(term, 'gene')) %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, HR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(term = gsub('gene', i, term),
           # add a column which states the model run
           model_covars = paste0(i,' +', i, ':', t, ' + ', t, ' + age_study_entry + sex + pc1-4 + stage_grouped + sact_over6weeks_before_studyentry + current_sact')
    )
  }
}
# merge into a single df and look at the p values
results_df <- lapply(results, function(x) do.call(rbind, x)) %>%
  do.call(rbind,.) %>%
  # just extract interaction terms
  filter(str_detect(term, ':')) %>%
  arrange(term, pval) %>%
  mutate(description = 'cancer_type_binary uses all_other_tumour_types as the reference')

results_df %>% filter(pval<0.05)
write.table(results_df, './results/interaction_terms_tumour_gene.tsv', 
            row.names=F, quote = F, sep = '\t')

# Time interactions ####
# check  Cox PH assumptions for the top FDR signifcant genes 
# for this I need to store the complete CoxPH model
fdr_sig_genes <- fread('./results/cox_genes_minadj.tsv') %>% filter(fdr_p <= 0.1) %>% pull(gene)
print(fdr_sig_genes)
# here I include gene:time_to_status in the model (linear interaction term between the gene and time)
# create vector of unique timepoints to split over:
cut1 <- 365*1
cut2 <- 365*2
cut3 <- 365*3
cut4 <- 365*4

res <- list()

for (i in fdr_sig_genes) {
  variants_by_gene <- all_variants %>% filter(gene == i)
  temp_data <- ac %>% 
       # identify whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  # split each record according to time_cuts. 
  newdata <- survSplit(Surv(time = time_to_status, event = binary_status) ~ 
                         gene + age_study_entry + sex + pc1 + pc2 + pc3 + pc4, data=temp_data, 
                       cut=c(cut1, cut2, cut3, cut4), id = 'participant_id') %>%
    # create a new variable time_cut
    mutate(time_point = as.factor(case_when(time_to_status <= cut1 ~ 1,
                                            time_to_status <= cut2 ~ 2,
                                            time_to_status <= cut3 ~ 3,
                                            time_to_status <= cut4 ~ 4,
                                            time_to_status > cut4 ~ 5)))
  # run the model with the linear time interaction term
  res[[i]] <- coxph(Surv(tstart, time_to_status, event = binary_status) ~ 
                      gene + time_point + time_point:gene + age_study_entry + 
                      sex + pc1 + pc2 + pc3 + pc4, data=newdata) %>%
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4 and disease_type)
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    dplyr::rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    dplyr::select(term, HR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(term = gsub('gene', i, term))
}


res_df <- do.call(rbind, res) %>% 
  # select the results relating to gene name, not the other covars
  filter(str_detect(term, str_c(fdr_sig_genes, collapse = '|'))) %>%
  # be explicit about which model was run
  mutate(model = 'gene + time_point + time_point:gene + age_study_entry + sex + pc1 + pc2 + pc3 + pc4')

write.table(res_df, './results/cox_genes_time_interaction.tsv', row.names = F, quote =F, sep = '\t')

## time stratified ####
# run stratified analyses for each time period and compare with interaction results
res <- list()

for (i in fdr_sig_genes) {
  variants_by_gene <- all_variants %>% filter(gene == i)
  temp_data <- ac %>% 
        # identify whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  # split each record according to time_cuts. 
  newdata <- survSplit(Surv(time = time_to_status, event = binary_status) ~ gene 
                       + age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                       data=temp_data, cut=c(cut1, cut2, cut3, cut4)) %>%
    # create a new variable time_cut
    mutate(time_point = as.factor(case_when(time_to_status <= cut1 ~ 1,
                                            time_to_status <= cut2 ~ 2,
                                            time_to_status <= cut3 ~ 3,
                                            time_to_status <= cut4 ~ 4,
                                            time_to_status > cut4 ~ 5)))
  for (t in 1:5) { 
    newdata_timepoint <- newdata %>% filter(time_point == t)
    # run the model with the time interaction term, with the data stratified by time point
    res[[i]][[t]] <- coxph(Surv(tstart, time_to_status, event = binary_status) ~ gene + age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                           data=newdata_timepoint) %>%
      tidy(exp=T, conf.int = T) %>% 
      # round estimates to 4 decimals
      mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
      # the output provides the HR for VTE associated with each of the exposure variables in the formula,
      # after adjusting for all the other variables
      # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4 and disease_type)
      # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
      dplyr::rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
      # select relevant columns for output
      dplyr::select(term, HR, l95ci, u95ci, pval) %>%
      # add a column which states the model run
      mutate(model_covars = paste(i, '+ age_study_entry + sex + pc1 + pc2 + pc3 + pc4')) %>%
      # replace 'gene' with the actual name of the gene queried
      mutate(term = gsub('gene', i, term)) %>%
      # add the timepoint being queried
      mutate(time_point = t) %>%
      # add number of VTE cases occurring in that time period
      mutate(ncase = newdata_timepoint %>% filter(binary_status == 1) %>% nrow)
  }
}

# merge into single df
res_df <- lapply(res, function(x) do.call(rbind, x)) %>%
  do.call(rbind,.) %>%
  # retain only the results relating to the genes
  filter(term %in% fdr_sig_genes) %>%
  # arrange in order of term and then timepoint
  arrange(term, time_point) %>%
  #reorder columns
  relocate(term, time_point, ncase)

#save
write.table(res_df, './results/cox_genes_time_stratified.tsv', row.names = F, quote =F, sep = '\t')
  
# Ancestry-stratification #####
  # first create vector of ancestry groups
  ancestries <- ac %>% 
    # create vector of genetically inferred ancestries (from PCs)
    pull(genetically_inferred_ancestry_nothr) %>% unique
  
  results <- list()
  # for each disease type run the analysis
  for (a in ancestries) {
    # filter the analysis cohort for ancestry group a
    ac_ancestry_stratified <- ac %>% filter(genetically_inferred_ancestry_nothr == a)
    # create a new list to store results for that ancestry
    results[[a]] <- list()
    # now loop through the gene list to run the following function on each gene
    for (i in fdr_sig_genes) {
      # filter the all_variants table for the gene in question
      variants_by_gene <- all_variants %>% filter(gene == i)
      # create a new column in the temp_data frame
      temp_data <- ac_ancestry_stratified  %>% 
        #  identifies whether the participants sample appears in the variants_by_gene table
        mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
      # now run the cox model 
      results[[a]][[i]] <- coxph(Surv(time = time_to_status, event = binary_status) ~ gene + 
                                   age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                                 data = temp_data) %>%
        # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
        tidy(exp=T, conf.int = T) %>% 
        # round estimates to 4 decimals
        mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
        # the output provides the HR for VTE associated with each of the exposure variables in the formula,
        # after adjusting for all the other variables
        # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4 and disease_type)
        filter(term == 'gene') %>% 
        # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
        rename('gene' = term, 'HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
        # select relevant columns for output
        select(gene, HR, l95ci, u95ci, pval) %>%
        # replace 'gene' with the actual name of the gene queried
        mutate(gene = i, ancestry = a)
    }}
  
  # note had issues with model convergence for some ancestry-gene combinations due to sparse data
  results_df <- lapply(results, function(x) do.call(rbind, x)) %>%
    do.call(rbind,.) %>%
    arrange(ancestry, pval) %>%
    # add a column which states the model run
    mutate(model_covars = 'gene + age_study_entry + sex + pc1-4') %>%
    # a few of the subgroups have uci95 = Inf as the model did not converge due to sparse data
    # for these rows the analysis is invalid therefore substitute with NA vals
    mutate(HR = case_when(u95ci == 'Inf' ~ NA, .default = HR),
           l95ci = case_when(u95ci == 'Inf' ~ NA, .default = l95ci),
           pval = case_when(u95ci == 'Inf' ~ NA, .default = pval),
           u95ci = case_when(u95ci == 'Inf' ~ NA, .default = u95ci)
    )
  
  
  write.table(results_df, './results/cox_genes_ancestry_stratified.tsv', row.names=F, quote = F, sep = '\t')
 
 
  # Anti-coagulation-stratification #####
  
  results <- list()
  # for each of prior anticoagulation indication == Y vs N run the analysis
  for (a in c('yes', 'no')) {
    # filter 
    ac_anticoag_subgroup  <- ac %>% filter(prior_anticoag_indication == a)
    # create a new list to store results for each group
    results[[a]] <- list()
    # now loop through the gene list to run the following function on each gene individually
    for (i in fdr_sig_genes) {
      # filter the all_variants table for the gene in question
      variants_by_gene <- all_variants %>% filter(gene == i)
      # create a new column in the temp_data frame
      temp_data <- ac_anticoag_subgroup  %>% 
        #  identifies whether the participants sample appears in the variants_by_gene table
        mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
      # now run the cox model 
      results[[a]][[i]] <- coxph(Surv(time = time_to_status, event = binary_status) ~ gene + 
                                   age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                                 data = temp_data) %>%
        # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
        tidy(exp=T, conf.int = T) %>% 
        # round estimates to 4 decimals
        mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
        # the output provides the HR for VTE associated with each of the exposure variables in the formula,
        # after adjusting for all the other variables
        # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4 and disease_type)
        filter(term == 'gene') %>% 
        # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
        rename('gene' = term, 'HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
        # select relevant columns for output
        select(gene, HR, l95ci, u95ci, pval) %>%
        # replace 'gene' with the actual name of the gene queried
        mutate(gene = i, prior_anticoag_indication = a)
    }}
  

  results_df <- lapply(results, function(x) do.call(rbind, x)) %>%
    do.call(rbind,.) %>%
    # add a column which states the model run
    mutate(model_covars = 'gene + age_study_entry + sex + pc1-4')
  
  write.table(results_df, './results/cox_genes_anticoag_stratified.tsv', row.names=F, quote = F, sep = '\t')
  
  
  # New time Zero ######
  # sensitivity analysis with time zero at cancer diagnosis
  # This cohort includes people who experienced a VTE between diagnosis and study entry
  # time zero has been reset at diagnosis (instead of study entry)
ac_sensitivity_diag_equals_time0 <- fread('./data/ac_sensitivity_diag_equals_time0.csv') %>%
    # create a binary status column for cox models where death is censored
    mutate(binary_status = case_when(status == 1 ~ 1, .default = 0),
    # calculate age study entry
    age_study_entry = as.integer(format(study_entry, "%Y")) - year_of_birth)
  
## variants time0 cohort ####
  #NB the combined_variants file was filtered for the original analysis cohort
  # so also have to remake this file
  snvs_qc_filtered <- fread('./data/snvs_vaf_5percent_lof_genesubset_qc_filtered.csv') %>% 
    select(tumour_sample_platekey, gene)
  large_SVs_filtered <- fread('./data/large_SVs_filtered.csv') %>% 
    select(tumour_sample_platekey, gene_query) %>%
    rename('gene' = 'gene_query')
  # combine the list of snvs and large SVs
  all_variants_sensitivity_time0 <-rbind(snvs_qc_filtered, large_SVs_filtered) %>% 
    # as only the gene name and tumour_sample_platekey are represented,
    # participants will only be counted once for each gene even if they carry multiple different variants
    distinct %>%
    # this time filter participants in the analysis_sensitivity cohort
    filter(tumour_sample_platekey %in% ac_sensitivity_diag_equals_time0$tumour_sample_platekey)
  rm(large_SVs_filtered, snvs_qc_filtered)
  

# retrict to new untreated cancer: pts entering study within 42 days of diagnosis who had not received prior treatment #####
df_new_untreated <- ac_sensitivity_diag_equals_time0 %>% 
  filter(new_untreated_cancer == 'yes')

# remove the old df so there is no confusion
rm(ac_sensitivity_diag_equals_time0)

count(df_new_untreated) %>% 
  write.table(.,
              './results/N_time0_new_untreated_sensitivity_analysis.txt', row.names = F)
#loop through the gene list to run the following function on each gene individually
results <- list()
for (i in gene_list) {
  # filter the all_variants table for the gene in question
  variants_by_gene <- all_variants_sensitivity_time0 %>% filter(gene == i)
  # create a new column in the analysis cohort dataframe 
  temp_data <- df_new_untreated %>%
    # which identifies whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  # now run the cox model
  results[[i]] <- coxph(Surv(time = time_from_diagnosis_to_status, event = binary_status) ~ gene + age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                        data = temp_data) %>% 
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
    filter(term == 'gene') %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    rename('gene' = term, 'HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(gene, HR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(gene = i)
}

# merge into a single df and look at the fdr adjusted p values
results_df_newuntreated <- do.call(rbind, results) %>%
  # add adjusted p values
  mutate(fdr_p = p.adjust(pval, method = 'fdr')) %>% 
  arrange(fdr_p) %>%
  # add a column which states the model run
  mutate(model_covars = paste(gene, '+ age_study_entry + sex + pc1-4'),
         sensitivity_analysis = 'new-untreated; time zero as cancer diagnosis')

fwrite(results_df_newuntreated, './results/cox_new_untreated_timezero_diag_genes.csv',
       row.names = F)

###full adj ######
df_new_untreated <- ac_sensitivity_diag_equals_time0 %>% 
  filter(new_untreated_cancer == 'yes')

tdata_sensitivity_diag_equals_time0 <- df_new_untreated %>% 
  # recalculate time from diagnosis to SACT (not time from study entry)
  mutate(time_diagnosis_to_sact = if_else(
    diagnosisdatebest < study_entry,
    # diagnosis before entry: add the delay
    time_diag_to_biopsy + time_to_current_sact,  
    # otherwise diagnosis after/at entry: keep original
    time_to_current_sact), 
    # create the stage_grouped column
    stage_grouped = case_when(stage_numeric_imputed <= 2 ~ 'early',
                              stage_numeric_imputed > 2 &
                                stage_numeric_imputed <= 4 ~ 'late',
                              is.na(stage_numeric_imputed) ~ 'unknown'),
    stage_grouped= factor(stage_grouped, levels = c('early', 'late', 'unknown')),
    # add +0.5 to time_to_status only where time_to_status =0
    # this step is required becuse tmerge does not accept stop times of 0
    # see vignette by Therneau et al
    time_from_diagnosis_to_status = case_when(time_from_diagnosis_to_status == 0 ~ 0.5,
                                              .default = time_from_diagnosis_to_status)
  )
# now perform tmerge to split f/u times for participants who received sact after study entry
sdata_sensitivity_diag_equals_time0 <- tmerge(
  data1 = tdata_sensitivity_diag_equals_time0,
  data2 = tdata_sensitivity_diag_equals_time0,
  vte = event(time_from_diagnosis_to_status, binary_status), 
  current_sact = tdc(time_diagnosis_to_sact), 
  id = participant_id,
  options= list(idname="participant_id")) %>%
  # select only relevant columns
  select(participant_id, tumour_sample_platekey, age_study_entry,
         sex, pc1, pc2, pc3, pc4,
         disease_type, stage_grouped,
         current_sact,
         tstart, tstop, vte)

results <- list()
#loop through the gene list to run the following function on each gene individually
for (i in gene_list) {
  # filter the all_variants table for the gene in question
  variants_by_gene <- all_variants_sensitivity_time0 %>% filter(gene == i)
  # create a new column in the analysis cohort dataframe NB using sdata here, 
  # where participants who receive sact after study entry have split time intervals
  temp_data <- sdata_sensitivity_diag_equals_time0 %>% 
    #  identifies whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  # now run the cox model
  results[[i]] <- coxph(Surv(tstart, tstop, vte) ~ gene + 
                          age_study_entry + sex + pc1 + pc2 + pc3 + pc4 +
                          disease_type + stage_grouped + 
                          # don't need to include SACT Before study entry because
                          # in this sensitivity cohort no one had SACT
                          current_sact, 
                        data = temp_data) %>%
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4 and disease_type)
    filter(term == 'gene') %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    rename('gene' = term, 'HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(gene, HR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(gene = i,
           model_covars = paste(i, '+ age_study_entry + sex + pc1-4 + disease_type + stage_grouped + current_sact')
    )
}

# merge into a single df and look at the fdr adjusted p values
results_df <- do.call(rbind, results) %>%
  arrange(pval) %>%
  # add a column which states the model run
  mutate(sensitivity_analysis = 'new-untreated; time zero as cancer diagnosis; fully_adj')


write.table(results_df, './results/cox_new_untreated_timezero_diag_genes_fulladj.csv',
            quote = F, row.names = F, sep = '\t')



# End ####