# Filter SNVs from cancer_tiered_variants_extract.sh and  filter large_structural_variants
# N Cornish

# Config #####
# First run the config.R script to set up libraries and the labkey_sql function
setwd(Sys.getenv('wd'))
source('./scripts/config421.R')

# SNVs ####
# these come from the cancer_tier_and_domain table and were extracted with cancer_tiered_variants_extract.R
# these have been downloaded for ALL participants. Contains tiered somatic variants in >19,000 canonical transcripts
# see the cancer_tiered_variants_extract.R for details
# read into local R studio
sample_info <- fread('./data/cancer_tiered_variants_samples.csv') %>%
  # filter only for somatic variants (note some germline variants appear in this table)
  filter(origin == 'somatic')
vep_info <- fread('./data/cancer_tiered_variants_vep.csv')


## COSMIC ####
# load COSMIC: /public_data_resources/COSMIC/v95/CosmicMutantExport.tsv.gz 
# filter this file in terminal (bash) to only include variants marked as PATHOGENIC (FATHMM prediction) in COSMIC
# the filtered file is here
cosmic <- fread('./data/CosmicMutantExport.tsv') %>% 
  # select only relevant columns (don't need COSMIC mutation ID etc as this parameter is not present in vep_info or sample_info)
  # to merge with vep_info need the GrCh38 coordinates and the chr, position and bp change (V9)
  select(V1, V2, V5, V6, V7, V8, V9) %>%
  # the colnames have not been preserved therefore add these back in
  setnames(., colnames(.), c('gene_id', 'transcript_id', 'cds_change', 
                             'GrCh', 'chr_pos', 'FATHMM_consequence', 'HGVS')) %>%
  # in order to merge with the vep info file it is necessary to first re-format the transcript_id and chr_pos so they align
  # the cosmic transcript ids are in the format ENTS****.$ 
  #(remove the digit after the decimal as the vep_info transcripts do not have this)
  mutate(transcript_id = gsub("\\.\\d+$", "", transcript_id))
                            

# before proceeding check that all coordinates are in GrCh38
table(cosmic$GrCh, useNA = 'always')
# check that all mutations in this file are marked as PATHOGENIC (if filtering has worked correctly this will be 'yes')
table(cosmic$FATHMM_consequence, useNA = 'always')

# Now extract just the HGVS notation along with the transcript for all the pathogenic Cosmic mutations
cosmic_pathogenic <- cosmic %>% 
  # concatenate the HGVS notation with the transcript ID
  mutate(transcript_HGVS = paste0(transcript_id, '_', HGVS)) %>%
           pull(transcript_HGVS) %>% unique


# Ensembl gene list (see the script gene_list for how this was produced)
ensembl_protein_genes <- fread('./data/ensembl_protein_genes_annotated.csv')

### sample_level ####
# first filter the cancer_tier_domain_variants ONLY for protein coding genes
sample_protein_coding <- sample_info %>% distinct %>% filter(gene %in% ensembl_protein_genes$gene_name)
# from sample_info, remove rows where the vaf is <5% (low tumour purity or low confidence in the call)
sample_vaf_5percent <- sample_protein_coding %>% 
  filter(vaf >= 0.05) 

### vep_level ####
# some variants are represented multiple times in the vep table
# this is due to the fact that the columns: cellbase_consequence, clinical_significane and relevance, sometimes have 2 slightly different values for the same variant / transcript
# from vep_info remove variants which are not flagged as potentially function altering in the relevance column as this combined info from cellbase_consequence and clinical_significance
vep_filtered <- vep_info %>%
  # create a new column with transcript_HGVS notation so that variants can be matched with cosmic_pathogenic
  # HGVS notation is as follows chr:g.posREF>ALT
  mutate(HGVS = paste0(transcript, '_', chr, ':', 'g.', pos, ref, '>', alt)) %>%
  # use this column to match with the cosmic_pathogenic list
  mutate(cosmic_pathogenic = case_when(HGVS %in% cosmic_pathogenic ~ T, .default = F)) %>%
    # remove the cellbase_consequence and clinical_significance columns and select unique entries
    # now filter only for variants which are marked as LOF on cellbase, 
    # OR pathogenic/likely pathogenic in ClinVAr OR pathogenic in COSMIC
  filter(relevance == '(likely)pathogenic'| # means represented as potentially pathogenic in clinvar
           relevance == 'LoF'| # in silico cellbase prediction
           relevance == 'path_LoF'| # both insilico cellbase prediction and in clinvar
           cosmic_pathogenic == T # pathogenic in COSMIC
         ) %>% 
  # since I have already filtered on relevance, 
  # select unique entries while ignoring the values in cellbase_consequence, clinical_significance and relevance 
  # in variants which are duplicated across rows the first row is retained
  distinct(select(., -c(relevance, cellbase_consequence, clinical_significance, cosmic_pathogenic)), 
           .keep_all=TRUE) # keep all columns


### merge filtered sample file with the filtered vep file ####
# Use inner join so that only variants which are represented in the cohort AND in the vep_filtered criteria are retained)
snvs_vaf_5percent_lof <- inner_join(sample_vaf_5percent, vep_filtered, 
                      by = c('chr', 'pos', 'ref', 'alt', 'transcript'))
 
### save #### Note this file contains SNVs all participants in cancer_analysis 
fwrite(snvs_vaf_5percent_lof, './data/snvs_vaf_5percent_lof.csv', row.names = F)

# this output file is then used as an input for flagged_snv_submit_v3.sh 
# in order to perform further filtering based on VCF quality flags

# clear environment
rm(sample_info, vep_info, sample_vaf_5percent, vep_filtered, 
   sample_protein_coding, cosmic, flow_diag, cellbase_lof, COSMIC_pathogenic, clinvar_pathogenic, cosmic_pathogenic)


# Large SVs #####
# these have been extracted using large_SV_query_v4.sh and formatted using /scripts/large_SV_format_v2.R
# in these scripts I used the getSVCNVperGene R package to pull all manta called somatic variants from the json files
# this was a chunked extraction due to file size. The foramtting script simply merges them all into a single file
# the ensembl gene list was used as a starting point

large_SVs <- fread('./data/large_SV_all_somatic.csv') %>% distinct

large_SVs_filter1 <- large_SVs %>%
  #  remove any SVs of size < 50; keep NA values as the translocations have NA for size
  filter(size >= 50|is.na(size))

## save ###
fwrite(large_SVs_filter1, './data/large_SVs_filtered.csv', row.names = F)


