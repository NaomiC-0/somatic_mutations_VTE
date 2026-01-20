# R script for extracting all variants from the cancer_tier_and_domain variants table on labkey
# N Cornish

# Query cancer-tiered variants table ####
# get a list of all the variants in the cancer_tiering table
# this contains counts of all participants with a small genetic variant of moderate or high impact
# variants present in < 5 participants are not reported
# description of how this report was produced are found here: https://re-docs.genomicsengland.co.uk/gene-centric-snv-report/
# this includes non-synonymous, splice site and RNA gene small variants for all participants
# cellbase_consequence: from Cellbase
# clinical_significance: annotation from ClinVar (version 2022-08-2024)
# relevance: LoF if the variant is marked in cellbase as deleterious (as per below); pathogenic if the variant is marked as pathogenic in clinvar and path_LOF if both true
# transcript ablation, splice acceptor variant, splice donor variant, stop gained, frameshift, stop lost, start lost, inframe insertion, inframe deletion 
# domain 1 = virtual panel of potentially actionable genes; 
# domain 2 = virtual panel of cancer-related genes curated by Sanger's Cancer Gene census
# domain 3: = any other genes

source('./scripts/config421.R')

table_queried <- "cancer_tier_and_domain_variants" 

# Attempting to load this whole table into R generally causes errors as the file is large
# Therefore do in 2 steps
# Step 1: download columns relating to the individual samples
columns_selected <- c("tumour_sample_platekey", "origin", "gene", "transcript", "chr", "pos", "ref", "alt", "vaf")
sample_info <- labkey_select()
# Step 2: download columns relating to the vep-annotation
columns_selected <- c("chr", "pos", "ref", "alt", "transcript", "cellbase_consequence", "clinical_significance", "relevance")
# ask for distinct rows for the vep info (avoid duplicate entries about the same transcript that have been associated with different participants)
vep_info <- labkey_select() %>% distinct()

# save as 2 separate files initially
fwrite(sample_info, './data/cancer_tiered_variants_samples.csv', row.names = F)
fwrite(vep_info, './data/cancer_tiered_variants_vep.csv', row.names = F)
# formatting of these files is shown in 'filter_somatic_variants.R' script 