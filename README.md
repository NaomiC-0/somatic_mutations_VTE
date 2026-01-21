---
title: "Workflow for: Associations of tumour somatic mutations with cancer-associated venous thromboembolism"
format: html
editor: visual
---

# Introduction

This code accompanies the manuscript titled:

'Associations of tumour somatic mutations with cancer-associated venous thromboembolism'.

DOI ....

All scripts have gone through the Genomics England (GEL) airlock approval process. Some scripts were adapted from training materials provided by GEL. See: <https://re-docs.genomicsengland.co.uk/tutorials/> and <https://re-docs.genomicsengland.co.uk/workflows/>

A config file `./config421.R` with file-paths to relevant directories, packages and bespoke functions for use in the GEL research environment (RE) is required. Scripts for pulling and formatting genetic data relevant to this analysis are specific to the GEL data structure and can also be found within the research environment.

Analyses were run using R version 4.2.1. Genetic data was queried using bedtools v2.31, bcftools v1.21 and plink v2.0 (see relevant scripts)

```{r}
sessionInfo()
```

```         
R version 4.2.1 (2022-06-23)
attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] getSVCNVperGene_2.0.0 RSQLite_2.3.7         forcats_1.0.0         survminer_0.4.9      
 [5] ggpubr_0.6.0          Rlabkey_3.2.2         jsonlite_1.8.9        httr_1.4.7           
 [9] kableExtra_1.4.0      muhaz_1.2.6.4         tibble_3.2.1          minty_0.0.1          
[13] readODS_2.3.0         lubridate_1.9.3       tidyr_1.3.1           labelled_2.13.0      
[17] gt_0.10.1             gtsummary_1.7.2       ggplot2_3.5.1         readr_2.1.5          
[21] stringr_1.5.1         broom_1.0.6           survival_3.3-1        dplyr_1.1.4          
[25] data.table_1.15.4    

loaded via a namespace (and not attached):
 [1] nlme_3.1-157         cmprsk_2.2-12        bit64_4.0.5          numDeriv_2016.8-1.1 
 [5] tools_4.2.1          backports_1.5.0      utf8_1.2.4           R6_2.5.1            
 [9] metafor_4.6-0        DBI_1.2.3            colorspace_2.1-0     withr_3.0.0         
[13] tidyselect_1.2.1     gridExtra_2.3        bit_4.0.5            curl_5.2.1          
[17] compiler_4.2.1       cli_3.6.3            xml2_1.3.6           scales_1.3.0        
[21] survMisc_0.5.6       systemfonts_1.1.0    digest_0.6.37        rmarkdown_2.29      
[25] svglite_2.1.3        pkgconfig_2.0.3      htmltools_0.5.8.1    fastmap_1.2.0       
[29] rlang_1.1.5          rstudioapi_0.16.0    generics_0.1.3       zoo_1.8-12          
[33] vroom_1.6.5          car_3.1-2            magrittr_2.0.3       metadat_1.2-0       
[37] patchwork_1.2.0      Matrix_1.6-0         Rcpp_1.0.12          munsell_0.5.1       
[41] fansi_1.0.6          abind_1.4-5          lifecycle_1.0.4      stringi_1.8.4       
[45] carData_3.0-5        mathjaxr_1.6-0       grid_4.2.1           blob_1.2.4          
[49] parallel_4.2.1       crayon_1.5.2         lattice_0.20-45      haven_2.5.4         
[53] splines_4.2.1        hms_1.1.3            knitr_1.49           pillar_1.9.0        
[57] ggsignif_0.6.4       glue_1.8.0           evaluate_1.0.3       broom.helpers_1.15.0
[61] vctrs_0.6.5          tzdb_0.4.0           gtable_0.3.5         purrr_1.0.2         
[65] km.ci_0.5-6          cachem_1.1.0         xfun_0.50            xtable_1.8-4        
[69] rstatix_0.7.2        viridisLite_0.4.2    memoise_2.0.1        KMsurv_0.1-5        
[73] timechange_0.3.0 
```

# 1. Format cancer cohort

```{r}
source('format_cancercohort.R')
# output file: './data/pcs_imd_stage_cancer_ncras_treatment.csv'
```

This script pulls the cohort of participants who were recruited to the 100,000 Genomes Project cancer programme and who have ongoing consent.

It then performs the following steps

1\) Identifies participants with 'Gold standard' paired tumour and germline WGS (i.e. both WGS samples have passed Genomics England internal quality checks)

2\) For participants with duplicate tumour samples submitted for WGS, selects the earliest submitted sample, or if multiple samples submitted on the same date prioritises one sample using the following criteria: a) biopsy from primary tumour site b) highest tumour purity, c) preservation method (Fresh frozen over formalin-fixed parrafin embedded) d) WGS quality metrics.

3\) Links GEL records with clinical information provided by the National cancer registry (NCRAS). For patients with more than 1 NCRAS record, the record where diagnosis date is closest to date of tumour biopsy for WGS is retained. If duplicate records still persist, the NCRAS record which most closely matches the clinical information submitted at recruitment is retained

4\) Flags patients who do not have any linked NCRAS record or where the diagnosis in the linked NCRAS record differs significantly from the diagnosis recorded at GEL recruitment (based on ICD10 WHO version 2019 level 2 mappings). Patients who were recruited using stored/retrospective biopsies (taken prior to recruitment opening) are also flagged for exclusion.

5\) Identifies the first date when a participant received surgery (within 6 weeks prior to study entry or following study entry), and adds first date of systemic anticancer therapy (SACT) where applicable

6\) Adds cancer staging information where available.

# 2. Format treatment data

```{r}
source('sact_data.R')
#output file: './data/sact_formatted_cancercohort.csv'
```

**sact_data.R script**: links data from the Systemic Anti-cancer therapy database (SACT) which collates mandatory reports on cancer treatments from hospital trusts (2014 onwards).

Identifies the first date when a patient received SACT (within 6 weeks prior to or any time following study entry) and flags patients who received SACT \> 6 weeks before study entry.

# 3. Format VTE diagnoses

```{r}
# input file: './data/code_lists/VTE_ICD10_algorithm.csv' (see supplementary table in manuscript)
source('VTE_diagnoses_cancer.R')
# output file './data/hes_combined_VTE_cancercohort.csv'
```

These scripts select the first recorded diagnosis of VTE from the hospital episode statistics including: admitted patient care records (hes_apc), A&E records (hes_ae and ecds) and outpatient records (hes_op).

Censor dates are listed here: \[https://re-docs.genomicsengland.co.uk/release19]

Since ecds records use SNOMED codes, the ICD10 codes were mapped to relevant SNOMED codes using the code lists provided by the GEMINI project: <https://github.com/GEMINI-multimorbidity.>

# 4. Format analysis cohort

a\) Loads the relevant tables created in steps 1-3 and merges into single df

b\) Creates status and time_to_status columns:

Although some linked electronic health records available in the version 19 data release run until 2024, the latest date for treatment records is August 2022. Therefore we censored participants at 5 years from study entry, or on 31/7/2022 (whichever came soonest) as this is the latest date when all linked health data can reasonably be assumed to be complete.

One limitation of time to event modelling using HES data to note: for inpatient VTE diagnosis, the date in hes_apc reflects the date of admission during which the VTE occurred. This might not always be the exact date of the VTE

c\) Flags patients who should be excluded from analysis cohort based on genetic QC thresholds, data discrepancies between NCRAS and GEL data, retrospective tumour biopsies, missing essential covariates or VTE prior to study entry.

d\) identifies sub-cohorts for sensitivity analyses including `new_untreated_cancer` column: participants who had tumour biopsied for WGS within 42 days of initial cancer diagnosis [and]{.underline} did not have any SACT before study entry; `prior_anticoag_indication` column flags participants with a medical indication for anticoagulation before study entry (see supplementary table for list of relevant ICD10 codes; to pull these from Labkey use script `anticoag_indications.R`

e\) Creates a table of values which can be used for **STROBE flow diagram** showing how analysis cohort was derived. See Figure 1 of manuscript.

f\) Formats final analysis cohort table. Ensure column classes are correct (factor, date etc) and derive new columns where required so that the table can be fed into various different statistical models with minimal additional formatting. Also group 'HEAD_NECK' tumours together and group other rare tumour types as 'OTHER'; create a new 'binary_status' column where death is treated as a censoring event (for cox models). Re-code participants who received SACT/ surgery within 6 weeks prior to study entry as exposed to these treatments at baseline and group cancer stages into `early` vs `advanced` or `unknown` if staging not available. Select only the columns required for the analysis and ensure that columns are coded as the correct class

```{r}
# INPUT files
#'./data/pcs_imd_stage_cancer_ncras_treatment.csv'
# './data/sact_formatted_cancercohort.csv'
# './data/hes_combined_VTE_cancercohort.csv'
# './data/hes_combined_anticoag_cancercohort.csv' (see anticoag_indications.R)

source('format_analysis_cohort.R')

# OUTPUT  files
#'./results/STROBE_diagram_vals.tsv' --> Figure 1 in manuscript
#'./data/analysis_cohort_formatted.rds'
```

# 5. Cohort demographics

a\) Creates summary of cohort demographics as per table 1 of manuscript

b\) calculates VTE rate over time: Supplementary Figure 1

c\) creates a list of tumour subtypes in each disease_group: see supplementary table 2

```{r}
# input files
# './data/analysis_cohort_formatted.rds'
source('cohort_demographics.R')
# output files
# "./results/cohortsummary_bystatus.tsv" --> Table 2 in main manuscript
#"./results/cohortsummary__continuous_vars.tsv" --> text in main manuscript
# './results/hazard_function_vte_5yr.png' --> supplementary figure 1
# './results/summary_tumour_types_airlock.csv' --> supplementary table 2
```

# 6. Somatic mutations

Small variants were extracted from the cancer_tier_and_domain_variants table as described here: <https://re-docs.genomicsengland.co.uk/cancer_tiering/>. Large variants (including Indels/translocations/inversion \> 50bp are extracted from GEL JSON files which can be queried using the GEL R package: `getSVCNVperGene` as described here: <https://re-docs.genomicsengland.co.uk/somatic_sv/>. The JSON files are created from the NSv4 delivery which used Manta v0.28; the raw output has been filtered by an internal GEL tiering pipeline which annotates variants likely to impact transcript function. In the scripts below I query all protein coding genes in ensembl (\>17,000 genes in output)

```{r}
source('cancer_tiered_variants_extract.R')
#output files from this are 
# cancer_tiered_variants_samples.csv' and 
#cancer_tiered_variants_vep.csv'

#./scripts/large_SV_query.R'  -> best to integrate this into a .sh script and run on HPC
# this script pulls the structural variants from these JSON files for all patients in cancer_analysis for all protein coding genes in ensembl. CNVs are not included
#This script outputs a list of results in 51 chunks (stored as Rdata files): ./data/large_SV_temp/
# convert the Rdata files to csv files and merge into a single dataframe with this script: 
source('large_SV_format.R')
# delete temp files after running this script
# final output file  is
# './data/large_SV_all_somatic.csv'
```

ClinVar and Cellbase annotations are in the cancer_tier_and_domain_variants table. COSMIC variant annotations obtained from COSMIC database

```{bash}
# Use Cosmic v95 file
# this version uses FATHMM scores to assign pathogenicity to variants.
# to check positions of specific columns
zcat CosmicMutantExport.tsv.gz | head -1 | tr '\t' '\n' | nl | egrep "Gene name|Accession Number|Mutation CDS|GENOMIC_MUTATION_ID|MUTATION_ID|FATHMM prediction|GRCh|Mutation genome position|HGVSG"

# To only select variants which are marked as 'pathogenic' (best to run this on the HPC)
zgrep 'PATHOGENIC' CosmicMutantExport.tsv.gz | awk -F'\t' '{print $1, $2, $17, $19, $20, $25, $26, $29, $39}' > $wd/data/CosmicMutantExport.tsv
```

Filter small and large somatic variants. Retain small variants if a) predicted LoF according to Cellbase or b) likely pathogenic/pathogenic according to ClinVar or COSMIC AND VAF \>5% in that sample AND none of the below quality flags in VCFs. Retain all large variants \>50bp which are present in the GEL tiered JSON files as these have already been through a filtering pipeline.

```{r}
# input files are 
# './data/cancer_tiered_variants_samples.csv' 
#'./data/cancer_tiered_variants_vep.csv'
#'./data/large_SV_all_somatic.csv'
# script to filter variants as described in supplementary methods
source('./scripts/filter_somatic_variants.R')
# then 
system('./flagged_snv_submit_v3.sh') 
#(run on HPC) to perform further filtering of small variants based on VCF quality flags:
# CommonGermlineVariants (gAF >1%) (G);
# CommonGnomADVariant: variants with population germline AG>1% in gnomAD;
# ReccurrentSomaticVariants (sAF >5%) (R);
# BCNoise10Indel: variants in regions of sequencing noise (N);
# PONnoise50SNV:
# SomaticFisherPhred < 50, indicating somatic SNV is systematic mapping/sequencing error;
# SimpleRepeat: overlapping simple repeats (SR); HomopolimerIndel: variants intersecting reference homopolymers (H).
# output files from this are
# ./data/snvs_vaf_5percent_lof_genesubset_qc_filtered.csv
# ./data/large_SVs_filtered.csv
```

Combine small and large variants into a single file and calculate mutation frequencies for each gene

```{r}
# ac sample size (needed to calculate % with mutations in each gene)
ac <- readRDS('./data/analysis_cohort_formatted.rds')
ac_samplesize <- nrow(ac) %>% as.numeric()
# SNVs
snvs_qc_filtered <-fread('./data/snvs_vaf_5percent_lof_genesubset_qc_filtered.csv') %>% 
  select(tumour_sample_platekey, gene)
# SVs 
large_SVs_filtered <- fread('./data/large_SVs_filtered.csv') %>% 
  select(tumour_sample_platekey, gene_query) %>%
  rename('gene' = 'gene_query')

# merge them
combined_variants <-rbind(snvs_qc_filtered, large_SVs_filtered) %>% 
  # as only the gene name and tumour_sample_platekey are represented,
  # participants will only be counted once for each gene even if they carry multiple different variants
  distinct %>%
  # retain only participants in the analysis cohort
  filter(tumour_sample_platekey %in% ac$tumour_sample_platekey)

freq <- combined_variants %>%
    # ensure each participant is counted only once for each gene
  select(tumour_sample_platekey, gene_query) %>% distinct %>% 
  # count number and prop of participants with mutations in each gene
  group_by(gene_query) %>% summarise(no_participants=n(),
                               proportion = round((no_participants/ac_samplesize), digits=4)) %>%
  # round proportions and arrange in order of freq
  arrange(desc(proportion))

write.table(freq_combined, './results/frequency_combined_SNV_SV_genelevel.tsv', row.names=F, sep = '\t', quote = F)
```

# 7. Gene list for analysis

Create a list of Ensembl protein-coding genes to be analysed based on the a-priori protocol criteria: 1) gene present in COSMIC census 2) previously reported in VTE GWAS 3) mutated in 5% of cohort

```{r}
# input files: 
# ensembl release-109
# COSMIC_CENSUS https://cancer.sanger.ac.uk/cosmic/ (downloaded Jan 2025)
#./data/vte_unique_genes_sept24.ods (manually compiled list from published literature: see manuscript for citations)
# ./results/frequency_combined_SNV_SV_genelevel.tsv
source('./scripts/gene_list.R')
# output
# ./data/gene_list_for_analysis.txt
```

# 8. Tumour mutational burden and mutational signatures

This information is found in the cancer_analysis table on labkey. For each tumour sample, TMB is calculated as the total number of non-synonymous small somatic variants divided by the total length of coding sequence (32.61 Mb). Genomics England computes mutational signatures using the R package nnls Details here: https://re-docs.genomicsengland.co.uk/cancer_analysis/#tumour-mutational-burden-tmb\
Mutation signatures correspond to version2 of the COSMIC signatures (March 2015):<https://cancer.sanger.ac.uk/signatures/signatures_v2/>\

```{r}
# see https://re-docs.genomicsengland.co.uk/cancer_analysis/#tumour-mutational-burden-tmb
# outputs: 
# './data/tmb_and_mutsig_ac.csv'
# /results/mutational_signatures_stats.tsv
```

# 9. Statistical analyses

## a) Cox models

These scripts run the main (primary) mininally-adjusted analyses, fully adjusted analyses and sensitivity analyses as described in the manuscript

```{r}
# Clinical #####
# inputs
# analysis cohort dataframe: see ./scripts/format_analysis_cohort.R
source('coxmodels_clinical.R')
# outputs
# './results/cox_clinical_univar.tsv' -> supplementary table 4

# Genes #######
# inputs
# analysis cohort dataframe. see ./scripts/format_analysis_cohort.R
# gene list. see ./scripts/gene_list.R
# combined variants file: see section under somatic mutations
source('coxmodels_genes.R')
# outputs
#cox_genes_minadj.tsv --> supp table 5 and figure 2
#cox_genes_fulladj.tsv --> supp table 5 and figure 2
# cox_genes_ancestry_stratified.tsv --> supp table 6
#cox_genes_anticoag_stratified.tsv --> supp table 6
#cox_genes_time_interaction.tsv --> supp table 6
#cox_genes_time_stratified.tsv --> supp table 6
#cox_genes_tumour_stratified.tsv --> supp table 6 and figure 3
#cox_new_untreated_timezero_diag_genes.csv --> supp table 6

# TMB #######
# inputs
# analysis cohort dataframe:  see ./scripts/format_analysis_cohort.R
# tmb/signatures file: see scripts/get_TMB.R
source('coxmodels_TMB.R')
# outputs
#cox_TMB_minadj.tsv --> table 2 and figure 4 
#cox_TMB_fulladj.tsv --> table 2 and figure 4
#cox_TMB_ancestry_stratified.tsv --> supp table 6
#cox_TMB_anticoag_stratified.tsv --> supp table 6
#cox_TMB_time_interaction.tsv --> supp table 6
#cox_TMB_time_stratified.tsv --> supp table 6
#cox_TMB_tumour_stratified.tsv --> supp table 6 and figure 4
#cox_new_untreated_timezero_diag_TMB.csv --> supp table 6
# cox_TMB_leave_one_cancer_out.tsv --> supplementary figure 5

#  mutational signatures #######
# inputs
# analysis cohort dataframe:  see ./scripts/format_analysis_cohort.R
# tmb/signatures file: see scripts/get_TMB.R
source('coxmodels_signatures.R')
# outputs
#cox_signatures_minadj.tsv --> table 2 and figure 4 
#cox_signature_fulladj.tsv --> table 2 and figure 4
#cox_signature_ancestry_stratified.tsv --> supp table 9
#cox_signature_anticoag_stratified.tsv --> supp table 9
#cox_signature_time_interaction.tsv --> supp table 9
#cox_signature_time_stratified.tsv --> supp table 9
#cox_signature_tumour_stratified.tsv --> supp table 9 and figure 4
#cox_new_untreated_timezero_diag_signatures.csv --> supp table 9
# cox_signature_continuous_nonzeros.tsv --> supplementary figure 9

# PRS and thrombophilia variants #####
# inputs
# analysis cohort dataframe:  see ./scripts/format_analysis_cohort.R
# germline genotypes for FVL and PTG20210A from AggV2: /klarin_SNP_germline_genotypes_all_cancer_analysis.txt
# PRS scores for cohort: see plink_calculate_PRS.sh
#./data/participant_PRS_scores_klarin.txt.sscore
source('coxmodels_FVL_PT_carriers.R')
source('coxmodels_PRS.R')
# ouputs
#cox_PRS_minadj.tsv -->  supplementary table 10
# cox_PRS_fulladj.tsv --> supplementary table 10
# cox_PRS_gene_interactions.tsv -->  supplementary table 10, figure 4
#cox_PRS_signature_interactions.tsv --> supplementary table 10
#cox_PRS_TMB_interactions.tsv --> supplementary table 10
# cox_germlineSNP_somatic_genes.tsv --> supplementary table 10, figure 4
```

## b) competing risk models

```{r}
# inputs
# analysis cohort dataframe. see ./scripts/format_analysis_cohort.R
# gene list. see ./scripts/gene_list.R
# combined variants file: see section under somatic mutations
# tmb/signatures file: see scripts/get_TMB.R
source('./scripts/cmprsk_regression.R')
# outputs
# FG_genes_fulladj.tsv
# FG_genes_minadj.tsv
# FG_signatures_fulladj.tsv
# FG_signatures_minadj.tsv
# FG_TMB_fulladj.tsv
# FG_TMB_minadj.tsv
# cuminc_plot_CDKN2A.tsv
#cuminc_plot_KRAS.tsv
#cuminc_plot_TP53.tsv
#cuminc_plot_PCDH15.tsv

```
