# Gene list for analysis
# N Cornish 

# Config ####
# Script which creates a list of genes to take forward into analysis
setwd(Sys.getenv('wd'))
source('./scripts/config421.R')

# Ensembl gene list #####
# Used Ensembl version 109 (2023)
library(rtracklayer)
gff <- import('/public_data_resources/ensembl-data/release-109/Homo_sapiens.GRCh38.109.chr.gff3.gz')
gff_df <- as.data.frame(gff)
gff_df_genes <- gff_df %>% 
  # filter for genes
  filter(!is.na(gene_id)) %>% 
  # exclude mitochondrial genes
  filter(seqnames != 'MT') %>% 
  # seqnames is the chr; rename
  mutate(chr = as.character(seqnames)) %>%
  # convert to numeric where X and Y are 23 and 24 respectively
  mutate(chr = as.numeric(ifelse(chr == 'X', 23,
                      ifelse(chr == 'Y', 24, chr)))) %>%
  # select chr cooridnates, position, gene name, gene_id and biotype
  select(chr, start, end, Name, biotype, gene_id)

# save
fwrite(gff_df_genes, './data/ensembl_109_genes.csv', row.names = F)

# Annotate ensembl list ####
# read in file
ensembl <- fread('./data/ensembl_109_genes.csv') %>% rename('gene_name' = 'Name')

## cancer census ####
# identify the oncogenes and TSGs using the cancer gene census (v101) from https://cancer.sanger.ac.uk/cosmic/ (downloaded Jan 2025)
cancer_census <- fread('./data/COSMIC_CENSUS_v101.csv') %>%
# ensure blank character values coded as NA
  mutate(across(c(1:ncol(.)) & where(is.character), ~na_if(.,"")))
  # format column names. Unfortunately this table does not include Ensembl Gene IDs
  rename('gene' = 'Gene Symbol', 
         'role_in_cancer' = 'Role in Cancer') %>%
  mutate(cancer_census = case_when(is.na(role_in_cancer) ~ 'other', .default = role_in_cancer)) %>%
  # select releavnt columns
  select(gene, cancer_census)


## vte genes #####
# this includes all germline VTE loci reported by Lindstrom et al 2019, Klarin et al 2019, Thibord et al 2022, Ghouse et al, 2023 and He et al, 2023
vte_genes <- read_ods('./data/vte_unique_genes_sept24.ods') 

# merge with the ensembl protein list
ensembl_gene_list <- ensembl %>% left_join(., cancer_census, by = c('gene_name' = 'gene')) %>%
  # add a column which identifies whether it is a known VTE locus
  mutate(vte_gwas_locus = ifelse(gene_name %in% vte_genes$Gene, T, F)) %>%
  # filter for protein coding genes only
  filter(biotype == 'protein_coding')

fwrite(ensembl_gene_list, './data/ensembl_protein_genes_annotated.csv', row.names=F)

# Gene list for analysis ######
# Reload annotated ensembl genes
ensembl_gene_list <- fread('./data/ensembl_protein_genes_annotated.csv') %>%
  # ensure no blank values (they should all be NA)
  mutate(across(c(1:ncol(.)) & where(is.character), ~na_if(.,"")))

# Load the combined variant frequencies for this cohort (see step 6 of readme doc)
combined_vars_5percent <- fread('./results/frequency_combined_SNV_SV_genelevel.tsv') %>%
  filter(proportion >= 0.05)

# now subset the ensembl list only for variants which are in the cancer census, VTE gene or variant freq > 5%
gene_list1 <- ensembl_gene_list %>% 
  filter(gene_name %in% combined_vars_5percent$gene|!is.na(cancer_census)|vte_gwas_locus == T) %>%
  pull(gene_name) %>% unique

# Now identify whether there are at least 10 VTE cases per category (wildtype vs mutated). This is a power-pruning step declared a-priori.
# Load genetic variants ####
all_variants <- fread('./data/combined_variants.csv')
# this table only includes people in the analysis cohort
# each participant is represented once per gene

# Load analysis cohort
ac <- readRDS('./data/analysis_cohort_formatted.rds')

# NB this takes a really long time- better to run as a job on the HPC
# create 2x2 contingency tables showing VTE freq for every gene in gene_list1
results <- list()
for (i in gene_list1) { # take one gene at a time
  # filter the all_variants table for the gene in question
  variants_by_gene  <- all_variants %>% filter(gene == i)
  # create a new column in the ac data 
  temp_data <- ac %>% 
    # which identifies whether the participants sample appears in the variants_by_gene table
    mutate(gene = case_when(tumour_sample_platekey %in% variants_by_gene$tumour_sample_platekey ~ 1, .default = 0))
  # create a frequency table identifying the freq of VTE events
  # in people with a mutated vs unmutated copy of the gene
  results[[i]] <- table(genotype = temp_data$gene, 
                        vte_status = temp_data$binary_status #,
                        #cancer_type = temp_data$disease_type
                        ) %>% 
    as.data.frame %>%
    # create a new column with the name of the gene and name of the cancer (in this case 'PANCANCER')
    mutate(gene = i, 
           cancer = 'PANCANCER')}

results_df <- do.call(rbind, results)

write.table(results_df, './results/contingency_tables_pancancer_genes.csv', row.names= F, quote = F, sep = '\t')

# Filter gene-list1 ####
contin_tables <- fread('./results/contingency_tables_pancancer_genes.csv')

# Now subset genelist1 to ONLY genes where the minimum observation count is 10 in any given box of the 2x2 table
final_gene_list <- contin_tables  %>%
  group_by(gene) %>%
# first exclude any genes with less than 4 rows (these are the ones where genotype = 0 only i.e. no samples with somatic mutations)
  filter(n()>2) %>%
  # for each gene slice the minimum row (minimum no of observations in any box in the 2x2 table)
  group_by(gene) %>% slice_min(Freq) %>% ungroup %>%
  # filter only for genes where Freq>=10
  filter(Freq >= 10)

#save
final_gene_list %>% distinct(gene) %>% write.table(., './data/gene_list_for_analysis.txt', quote=F, row.names = F, sep = '\t')

# End of Script ####

