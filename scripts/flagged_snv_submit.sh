# Bash script for extracting the cancer tier and domain variants from the somAgg vcf on the hpc
# N Cornish. 
# adapted from the somAgg codebook: https://re-docs.genomicsengland.co.uk/somAgg/#code-book
# this .sh script has been streamlined from a longer Rmd workflow which I used in the RE but has not been retested;
# there may be some bugs! 
# It would be preferrable to put this into a pipeline (e.g. next-flow) if possible (breaking each step decsribed here into a new task)

#!/bin/bash  
#BSUB -q medium  
#BSUB -P <enter_project_code>  
#BSUB -o <file_path>/job.%J.out  
#BSUB -e <file_path>/job.%J.err  
#BSUB -J extract_tiered_snvs  
#BSUB -R "rusage[mem=100] span[hosts=1]"  
#BSUB -M 2  
#BSUB -n 1  
#BSUB -cwd <file_path>/data/flagged_snvs  

# CONFIG
export my_analysis_folder="${wd}/data/flagged_snvs/"
# list of pre-filtered somatic variants from cancer_tier_and_domain_variants
export my_variants="${path_to_my_variants_file}"
export final_output="${wd}/data/my_variants_qc_filtered.csv"
export variants="${my_analysis_folder}/my_variants.bed"
export path_to_chunknames="${path_to_somAgg_VCF}/additional_data/chunk_names/somAgg_chunk_names.bed"

# load modules
module load R/4.2.1
module load bedtools 
module load bcftools

# Step 1: Create .bed file of variants in R
R --vanilla
source('config421.R')
wd <- Sys.getenv('my_analysis_folder')
my_variants <- fread(Sys.getenv('my_variants'))
coords <- my_variants %>% # format chr columns to align with vcf format 1 -> 'chr1' etc
  mutate(chr = paste0('chr',chr),
         start = pos, end = pos) %>% 
  # select just variant coordinates
  select(chr, start, end) 
# write to a .bed file with no header  
write.table(coords, file=paste0(wd,'/my_variants.bed'), quote=F, row.names=F, col.names=F, sep = "\t")
# quit R 
q()

# Step 2: Use bed tools to locate correct vcf chunks
bedtools intersect -wo -a $variants -b $path_to_chunknames > ${my_analysis_folder}/variant_chunks.txt

# Step 3: Use output file from step 2 to derive vcf filepaths for each variants in my_variants list
R --vanilla
source('config421.R')
wd <- Sys.getenv('my_analysis_folder')
my_variants <- fread(Sys.getenv('my_variants')) %>% 
  setDT %>%  # convert to data.table
  # change the format of the chr column to match the chunks file
  mutate(chr = paste0('chr',chr)) %>% 
  # create a coordinates column in the format chr:pos
  mutate(coords = paste0(chr, ':', pos)) %>%
  # select only essential columns
  select(tumour_sample_platekey, coords, ref, alt) %>% distinct

chunks <- fread(paste0(wd, '/variant_chunks.txt')) %>% 
  # select just the chr, pos and filepath columns
  select(V1, V2, V9) %>% 
  setnames(., colnames(.), c('chr', 'pos', 'filepath')) %>%
  # create a coordinates column in the format chr:pos
  mutate(coords = paste0(chr, ':', pos)) %>% select(coords, filepath) %>%
  distinct %>% 
  # convert to data.table for efficient merging 
  setDT

# Use data.table to left join my_variants with chunks
my_variants_chunk_paths <- merge(my_variants, chunks, by = 'coords',
                                 # left join 
                                 all.x=T, all.y=F) %>%
  # place columns in following order: sample, coords, filepath
  relocate(tumour_sample_platekey, coords, filepath) %>% arrange(filepath)

# save with no header
fwrite(my_variants_chunk_paths, paste0(wd, '/variant_filepaths.csv'), col.names = F)

# quit R
q()

# STEP 4, use bcf tools to pull each of the my_variants from the relevant vcf chunk in a loop 
# NOTE this part is the slow part. It is not possible to run a lot of variants
# if variant_filepaths.csv has many rows it is better to run in chunks of up to 50-100 rows at a time
# this can be done by editting 'start' and 'end' for each chunk below
start=1
end=$(wc -l < variant_filepaths.csv)
# ensure working from ./data/flagged_snvs
for ((row = start; row <= end; row++)); do
nrow=$(wc -l < variant_filepaths.csv)
for ((row = 1; row <= nrow; row++)); do
    export sample=$(awk -F',' -v r="$row" 'NR==r {print $1}' variant_filepaths.csv)
    export variant=$(awk -F',' -v r="$row" 'NR==r {print $2}' variant_filepaths.csv)
    export vcf_filepath=$(awk -F',' -v r="$row" 'NR==r {print $3}' variant_filepaths.csv)
    export results="variant_flags.txt"

bcftools query -s $sample -r $variant -f '[%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/FILTER\t%SAMPLE\n]' $vcf_filepath >> $results

done

# Now filter the variants in R to identify which variants had QC flags in the vcf files
R
.libPaths(c(.libPaths(), "/tools/aws-workspace-apps/ce/R/4.2.1"))
library(dplyr)
library(data.table)

# set working directory variable
wd <- Sys.getenv('my_analysis_folder')
# read in the original variant file for comparison (to check have got all variants from the vcfs)
my_variants <- fread(Sys.getenv('my_variants')) %>% 
  # create a row_id comprised of variant coordinates (in the format tumour_sample_platekey_chrN:pos) 
  mutate(row_id = paste0(tumour_sample_platekey, '_', 'chr', chr, ':', pos))
# note there may be some duplicated row_ids where the same variant has been called across multiple transcripts for the same sample. 
#The row ID will simply be used to check that the variant has been correctly pulled from the somAgg file

# read in the variants which have been pulled from SomAgg
file_names <- list.files(path = wd, pattern = '*variant_flags*')
# if you have run the bcftools section in chunks this will read in all files and merge them into a single df (providing all filenames contain the string *variant_flags*)
somAgg_vars <- lapply(file_names, fread) %>% purrr::reduce(full_join) %>% 
  setnames(., colnames(.), c('chr', 'pos', 'ref', 'alt', 'FILTER', 'INFO', 'tumour_sample_platekey')) %>% 
  # create a unique_id column which matches my_variants
  mutate(row_id = paste0(tumour_sample_platekey, '_', chr, ':', pos)) %>%
  # if there is a  '.' in the FILTER or INFO column, replace with NA vals
  mutate(across(c(1:ncol(.)) & where(is.character), ~na_if(.,"."))) 

# check that all variants from my_variants have been correctly pulled from the somAgg VCFs
vars_missing <- my_variants %>% filter(!row_id %in% somAgg_vars$row_id)
if (count(vars_missing) !=0) {
  print(paste(as.numeric(nrow(vars_missing)), 'variants from the "my_variants" file have not been pulled from the somAgg vcfs:'))
  print(vars_missing %>% select(tumour_sample_platekey, chr, pos))} else {
    print('all variants from the "my_variants" file have been pulled from the somAgg vcfs')
  }
# if any variants are missing from the somAgg_vars need to go back and troubleshoot

# assuming count(vars_missing)=zero....continue
# note the low Q score and BCNoise flags are in the FILTER column 
# the other relevant flags (Common Germline variant etc) are in the INFO column
# merge the FILTER and INFO column from somAgg_vars with the my_variants file based on chr, pos, ref, alt and tumour_sample_platekey
# first ensure that chr and pos are coded the same for both files
somAgg_vars_reformat <- somAgg_vars %>%
  # remove the string 'chr' from somAgg_vars
  mutate(chr = gsub('chr', '', chr),
         # change pos to a character variable so class is the same as my_variants
         pos = as.character(pos)) %>%
  # exclude any duplicated rows
  distinct

# now left join with my_variants (so all variants in my_variants are retained)
my_variants_with_qc_info <- left_join(my_variants, somAgg_vars_reformat, by = c('row_id', 'chr', 'pos', 'ref', 'alt', 'tumour_sample_platekey'))

# Now filter the my_variants_with_qc_info according the following QC criteria as per Sosinsky et al: 
#a) remove vairants with germline allele freq > 1% in Genomics England or gnomAD; 
#b) remove recurrent somatic variants with a freq >5% in Genomics England dataset; 
#c) remove variants overlapping simple repears as defined by Tandem Repeats Finder; 
#d) remove small indels in regions with high levels of sequencing noise (BCNoise10Indel);
#e) remove variants with Fisher's Phred score < 50.
my_variants_qc_filtered <- my_variants_with_qc_info %>% 
  filter(FILTER != 'BCNoise10Indel' & 
           FILTER != 'LowQscore ' &
           FILTER != 'PONnoise50SNV' &
           FILTER != 'QSI_ref,BCNoise10Indel' &
  # note I have deliberately used != sign here to allow both FILTER == 'PASS' and NA values through 
  # also want to retain the is.na(INFO) - all other values in INFO are exclusions
           is.na(INFO)
           ) 

# call the output filepath which is set in config
output <- Sys.getenv('final_output') 
# save to my_variants_qc_filtered to the output filepath
fwrite(my_variants_qc_filtered, file = output, row.names = FALSE)

q()

# END OF SCRIPT
exit