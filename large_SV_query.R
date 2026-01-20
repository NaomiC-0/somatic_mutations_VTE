# Structural variants query
# Script for pulling large structural variants and CNVs
# Adapted from here: https://re-docs.genomicsengland.co.uk/somatic_sv/
# N Cornish
# this script takes a long time to run: submit it as a job on the HPC

# note the output after running this script for every protein coding gene in ensembl and then merging the chunks together is found is ./data/large_SV_all_somatic.csv
# see large_SV_format_v2.R for details
# delete temporary files afterwards to conserve storage space on RE

# config ####
setwd(Sys.getenv('wd'))
source('./scripts/config421.R')
library(RSQLite)
.libPaths(c(.libPaths(), "/tools/aws-workspace-apps/ce/R/4.3.3"))
library(getSVCNVperGene)

# this function uses a slightly different format for the database version to the config script
release_version = paste0('/',version)

# Get list of all samples from cancer_analysis table ####
table_queried <- "cancer_analysis"
columns_selected <- c("tumour_sample_platekey")
samples <- labkey_select() %>% unlist
length(samples)

# gene list ####
ensembl_list <- fread('./data/ensembl_protein_genes_annotated.csv')
ensembl_protein_genes <- ensembl_list %>% select(gene_name) %>% unlist

# split the ensembl list into chunks
ensembl_chunks <- split(ensembl_protein_genes, ceiling(seq_along(ensembl_protein_genes)/400))

for (sublist in seq_along(ensembl_chunks)) {
  genes_to_query <- ensembl_chunks[[sublist]]
# get SVs ####
#getSV queries structural variants.
# SVTYPE: BND= translocation; DEL=deletion, DUP=duplication; INV=inversion, INS=insertion (as defined by Manta)
# variant_domain = domain assiged to the variants by the tiering pipeline (domain1, 2, 3)
# split the ensembl gene list into chunks of 10 genes
gene_chunks <- split(genes_to_query, ceiling(seq_along(genes_to_query)/10))

# create empty list to store results; 
results <- vector("list", length(genes_to_query))  # fix the size of the list as the number of genes being queried
# Explicitly assign relevant gene names to the list before populating it with the loop results
# I added this step because with a long loop I found each element of the results list was named using a numeric operator
# in a smaller loop the gene name was retained even though the code was identical. 
names(results) <- genes_to_query

# loop
# for each chunk
for (chunk in gene_chunks) {
  # create a list to add results to
  chunk_res <- list()
  # for each gene in the chunk
  for(i in chunk) {
    # try running the Structural variant query
    chunk_res[[i]] <- try({
      getSV(gene = i, fusionOnly = F, diseaseType = NULL, participantID = NULL,
            # submit either list of participant IDs or platekeys, not both
            # version is the data version (set in config script)
            release_version = release_version, plateKey = samples)
    }, 
    # silent = T means error messages do not print to console
    silent = T)
    # if there is an error print this to the results list
    if (inherits(chunk_res[[i]], "try-error")) {cat("Error encountered for gene:", i, "\n")
      # then proceed with the next gene
      next
    }
  }
  # after running the chunk, map the output to the correct gene in the results list
  results[names(chunk_res)] <- chunk_res
  }

# save the list in a directory called large_SV_temp before proceeding with the next element in the ensembl list
results_name <- paste0('chunk_',sublist)
save(results, file=paste0('./data/large_SV_temp/large_SV_query_allensembl',results_name,'.Rdata'))
    }

# End of script
