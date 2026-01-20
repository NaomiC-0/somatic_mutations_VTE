# Script looking at interactions between Factor V Leiden / Prothrombin G20210A carrier status
# and somatic mutations
# sensitivity analysis
#N Cornish

rm(list=ls())
setwd(Sys.getenv('wd'))
source('./scripts/config421.R')
library(survival)
library(survminer)
# Load formatted analysis cohort ####
ac <- readRDS('./data/analysis_cohort_formatted.rds')

# Load germline genotypes ####
# read in the germline genotypes for all the SNPs in the Klarin PRS
# these were extracted and formatted using PRS_calculation_v3.Rmd
germline_genotypes <- fread('./results/klarin_SNP_germline_genotypes_all_cancer_analysis.txt') %>%
  # this is for whole cancer cohort: restrict just to analysis cohort
   filter(germline_sample_platekey %in% ac$germline_sample_platekey)

# Load somatic mutations ####
all_variants <-fread('./data/combined_variants.csv') %>% 
  # participants will only be counted once for each gene even if they carry multiple different variant
  distinct
fdr_sig_genes <- fread('./results/cox_genes_minadj.tsv') %>% filter(fdr_p <= 0.1) %>% pull(gene)

# additive model with each gene and PRS
freqs <- list()
results_SNP_only <- list()
results_SNP_mut_interac <- list()
for (s in c('rs6025', 'rs1799963'))  {
  # filter the genotype table just for that SNP
  genotype <- germline_genotypes %>% filter(rsid == s)
  # make a table of genotype frequencies. Note I already filtered just for ac above
  freqs[[s]] <- table(genotype$sample_genotype) %>% as.data.frame %>%
    dplyr::rename(genotype = Var1) %>%
    mutate(SNP = s) %>% relocate(SNP)
  # Now merge with the analysis cohort df
  tdata <- ac %>% 
    # merge with the germline genotype data
    left_join(., genotype, by = 'germline_sample_platekey') %>%
    # remove the people who have NA for genotype (just the ones in AggV2)
    filter(!is.na(sample_genotype)) %>%
    #code sample genotype as a factor
    mutate(sample_genotype = factor(sample_genotype),
           # add +0.5 to time_to_status only where time_to_status =0
           # this step is required becuse tmerge does not accept stop times of 0
           # see vignette by Therneau et al
           time_to_status = case_when(time_to_status == 0 ~ 0.5,
                                      .default = time_to_status ))
  # first look at associations between SNP carrier status and VTE
  # then run cox model
  results_SNP_only[[s]] <- coxph(Surv(time_to_status, binary_status) ~ sample_genotype + age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                   data = tdata) %>% 
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
    mutate(model_covars = paste0(s, '_genotype + age_study_entry + sex + pc1 + pc2 + pc3 + pc4'),
           # annotate the term with SNP ID so it is clear which SNP the result is for
           term = gsub('sample_genotype1', paste0(s,'_Het'), term),
           term = gsub('sample_genotype2', paste0(s,'_Hom'), term)) %>%
             # remove the other covariate terms as not relevant
             filter(str_detect(term, s))
    # Now look at SNP:gene interactions
  for (i in fdr_sig_genes) {
  # filter the all_variants table for the gene in question
    variants_by_gene <- all_variants %>% filter(gene == i)
  # create a new column in the analysis cohort dataframe 
  temp_data <- tdata %>% 
        # identify whether the participants sample appears in the variants_by_gene table
    mutate(gene = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0))
  # now run the cox model
  results_SNP_mut_interac[[s]][[i]] <- coxph(Surv(time = time_to_status, event = binary_status) ~ gene + sample_genotype + gene:sample_genotype +
                          age_study_entry + sex + pc1 + pc2 + pc3 + pc4, 
                        data = temp_data) %>% 
    # use broom package to extract coefficients, exponentiate to obtain HR and 95% CIs
    tidy(exp=T, conf.int = T) %>% 
    # round estimates to 4 decimals
    mutate_at(vars(estimate, conf.low, conf.high), ~round(., 4)) %>%
    # the output provides the HR for VTE associated with each of the exposure variables in the formula,
    # after adjusting for all the other variables
    # extract just the HR for the gene (after adjusting for age, sex and pcs 1-4)
    filter(str_detect(term, 'gene')|str_detect(term, 'genotype')) %>% 
    # rename columns for clarity. Note default CIs correspond to 95% CI as per package notes
    rename('HR' = estimate, 'l95ci' = conf.low, 'u95ci' = conf.high, 'pval' = p.value) %>%
    # select relevant columns for output
    select(term, HR, l95ci, u95ci, pval) %>%
    # replace 'gene' with the actual name of the gene queried
    mutate(term = gsub('gene', i, term))  %>%
    # replace 'sample_genotype' with the actual SNP
    mutate(term = gsub('sample', s, term)) %>%
    mutate(model_covars = paste0(i, ' + ', s, '+ ', i,':', s, ' + age_study_entry + sex + pc1-4'))
    }
  }

# merge into a single df and look at the fdr adjusted p values
results_SNP_df <- do.call(rbind, results_SNP_only)
results_SNP_mut_interac_df_rs6025 <- results_SNP_mut_interac[['rs6025']] %>% do.call(rbind,.)
results_SNP_mut_interac_df_rs1799963 <- results_SNP_mut_interac[['rs1799963']] %>% do.call(rbind,.)
results_SNP_mut_interac_df <- rbind(results_SNP_mut_interac_df_rs6025, results_SNP_mut_interac_df_rs1799963) %>%
  # remove the results for FVL homozygotes as low numbers and model rarely converges (unstable regression)
  filter(!str_detect(term, 'genotype2'))


rm(results_SNP_mut_interac_df_rs6025,results_SNP_mut_interac_df_rs1799963)

freq_df <- lapply(freqs, as.data.frame) %>% do.call(rbind, .) %>%
  # add the denominator and calculate 
  group_by(SNP) %>%
  mutate(prob = Freq / sum(Freq)) %>%
  ungroup()
# merge
results <- rbind(results_SNP_df, results_SNP_mut_interac_df)

write.table(freq_df, './results/SNP_freqs_FVL_PTG20210A.tsv', quote = F, row.names = F, sep = '\t')
write.table(results, './results/cox_germlineSNP_somatic_genes.tsv', quote = F, row.names = F, sep = '\t')

# PCDH15:FVL cuminc ####
# Explore interaction between PCDH15 and FVL

genotype <- germline_genotypes %>% filter(rsid == 'rs6025')
variants_by_gene <- all_variants %>% filter(gene == 'PCDH15')

tdata <- ac %>% 
  # merge with the germline genotype data
  left_join(., genotype, by = 'germline_sample_platekey') %>%
  # remove the people who have NA for genotype (just the ones in AggV2)
  filter(!is.na(sample_genotype)) %>%
  # for the purpose of this analysis exclude FVL homozygotes 
  # this is because there are very few homozygotes and we only want 2 groups for clarity in the cuminc plot
  mutate(FVL_status = case_when(sample_genotype == 2 ~ NA, .default = sample_genotype),
         # code PCHD15 status
           # identify whether the participants sample appears in the variants_by_gene table
         PCDH15_status = ifelse(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey, 1, 0),
         # create a combined category
        combined_category = case_when(
          FVL_status == 1 & PCDH15_status == 1 ~ 'FVL carrier and PCDH15 mutated',
          FVL_status == 1 & PCDH15_status == 0 ~ 'FVL carrier and PCDH15 wildtype',
          FVL_status == 0 & PCDH15_status == 1 ~ 'FVL wildtype and PCDH15 mutated',
          FVL_status == 0 & PCDH15_status == 0 ~ 'FVL wildtype and PCDH15 wildtype',
          .default = NA),
         # add +0.5 to time_to_status only where time_to_status =0
         # this step is required becuse tmerge does not accept stop times of 0
         # see vignette by Therneau et al
         time_to_status = case_when(time_to_status == 0 ~ 0.5,
                                    .default = time_to_status ))
# 
fg_data <- finegray(Surv(time_to_status_months, status) ~ .,
                    data = tdata, 
                    etype = 1,
                    id = participant_id)
# run the model. combined category is the PRS/PCDH15 combined category
model <- coxph(Surv(fgstart, fgstop, fgstatus) ~combined_category,
               data = fg_data)

# create the new_df for surv-fit
new_df <- data.frame(combined_category = as.factor(c('FVL carrier and PCDH15 mutated',
                                                     'FVL carrier and PCDH15 wildtype',
                                                     'FVL wildtype and PCDH15 mutated',
                                                     'FVL wildtype and PCDH15 wildtype'
)))
fit <- survfit(model, newdata = new_df)
# invert the fit to get cuminc instead of survival
fit$surv=1-fit$surv
# for more flexible altering
surv_data <- surv_summary(fit)
# map the correct term to each strata
levels(surv_data$strata) <- c('FVL carrier and PCDH15 mutated',
                              'FVL carrier and PCDH15 wildtype',
                              'FVL wildtype and PCDH15 mutated',
                              'FVL wildtype and PCDH15 wildtype')

write.table(surv_data, './results/cuminc_plot_FVL_PCDH15_cmprsk.tsv',
            row.names = F, quote = F, sep = '\t')

