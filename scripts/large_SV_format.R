# Structural variants formatting and filtering
# Naomi Cornish
# 06/03/2025

# this script takes output from large_SV_query.R and formats into a single df
setwd(Sys.getenv('wd'))
source('./scripts/config421.R')

# gene list ####
# reload the ensembl gene list that I used in the large_SV_query script
ensembl_protein_genes <- fread('./data/ensembl_protein_genes_annotated.csv') %>% 
  select(gene_name) %>% unlist
# split the ensembl list back into chunks
ensembl_chunks <- split(ensembl_protein_genes, ceiling(seq_along(ensembl_protein_genes)/400))

# get a list of all the chunked outputs from the large_SV_query_v3.R
setwd(Sys.getenv('wd'))
setwd('./data/large_SV_temp/files_renamed')
files <- list.files(pattern = 'Rdata')
# check all files are named as expected by extracting just the digits
files %>% gsub("\\D", "", .) %>% as.numeric %>% sort %>% print
# Convert results to df ####
# write any output from the console (e.g. error messages) to file
# ensure there is a directory called ./data/large_SV_temp/files_renamed/temp_dataframes/
sink(file = './temp_dataframes/console_output.txt')
# loop through each chunk in the files list (51 files corresponding to 51 ensembl chunks)
for (i in c(1:51)) {
  paste('running chunk', i) %>% print
  load(files[[i]]) # note the resulting environment variable is called 'results'
  # double check all genes in the ensembl chunk appear in the results
  genes_queried <- names(results)
  missing_genes <- ensembl_chunks[[i]][ensembl_chunks[[i]] %!in% genes_queried] 
  if(length(missing_genes > 0)){print('gene missing from output:')
    print(as.data.frame(missing_genes), row.names=F)}
  # check for NULL elements
  null_elements <- which(sapply(results, is.null))
  # if null elements are present write these to the console
  if(length(null_elements > 0)) {print('null elements present:')
    ensembl_chunks[[i]][null_elements] %>% as.data.frame %>% print(row.names=F)}
  # these sometimes occur due to genes which appear in the ensembl file with an ensembl_id but no gene name
  # also some NULL results if no participants with SVs in the gene
  
  # Remove any NULL elements from the results
  results_clean  <- results %>% Filter(Negate(is.null), .) %>%
    # results[!sapply(results, inherits, "try-error")] %>%
    # only keep genes where there is a result (i.e. where nrow >= 1)
    # so if there were no structural variants in the sample then do not write to table
    Filter(function(x) nrow(x) >=1, .)
  
  # add a new column to each element in the list 
  #which corresponds to the name of that element (i.e the gene queried)
  lapply(names(results_clean), function(name) {
    df <- results_clean[[name]]
    df$gene_query <- name
    df
  }) %>%
    # convert to a df
    do.call(rbind, .) %>%  
    # select the relevant columns
    select(participant_id, tumour_sample_platekey, variant_origin, gene_query, SVTYPE, 
           chromosome, start, end, size, gene_name_bp1, gene_name_bp2, ensembl_id, bp1_location,
           bp2_location, fusion_inferred) %>%
    # filter for somatic variants
    filter(variant_origin == 'somatic') %>%
    # save to file
    fwrite(., paste0('./temp_dataframes/chunk_',i,'.csv'), row.names = F)}
sink() # stop writing output to file

# Merge Dataframes ####
# first clear environment
rm(list=ls())
setwd('./temp_dataframes')
data_frame_names <- list.files(pattern = "*.csv")
# read in all the dataframes
results <- lapply(data_frame_names, fread)

# navigate back to wd
setwd(Sys.getenv('wd'))
# merge into a single table
do.call(rbind, results) %>%
  # save in the many data directory
  fwrite(., './data/large_SV_all_somatic.csv', row.names=F)

# delete anything in temp_dataframes to free up memory!



