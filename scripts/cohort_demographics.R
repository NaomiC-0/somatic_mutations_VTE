# Cohort demographics
# N Cornish

# Notes ######
# Script migrated from main workflow in order to streamline it for github
# It was originally part of a longer Rmd document and I have not re-run in this format
# Therefore there may be some bugs

ac <- readRDS('./data/analysis_cohort_formatted.rds')

# Baseline demographics #####
#total sample size
total_n <- nrow(ac) %>% as.numeric
n_no_vte <- ac %>% filter(binary_status == 0) %>% nrow %>% as.numeric
n_vte <- ac %>% filter(binary_status == 1) %>% nrow %>% as.numeric

# create a factor of the disease names so that in the summary table the disease_types appear in order of overall frequency
disease_orders <- ac %>% group_by(disease_type) %>% summarise(count = n()) %>% arrange(desc(count)) %>% pull(disease_type) %>% factor(., levels=unique(.))

# create a temp df where disease_type is arranged in disease_orders above and use this in the matrix
temp <- ac %>% mutate(disease_type = factor(disease_type, levels = disease_orders))

matrix <- rbind(
  'sex' = rep(NA, 2),
  table(ac$sex, ac$binary_status), 
  'age_category' = rep(NA, 2),
  table(ac$age_category_studyentry, ac$binary_status),
  'ancestry' = rep(NA, 2),
  table(ac$genetically_inferred_ancestry_nothr, ac$binary_status),
  'cancer' = rep(NA, 2),
  table(temp$disease_type, temp$binary_status),
  'stage of tumour' = rep(NA, 2),
  table(ac$stage_numeric_imputed, ac$binary_status),
  'stage category' = rep(NA, 2),
  table(ac$stage_grouped, ac$binary_status),
  'tumour sample' = rep(NA, 2),
  table(ac$tumour_type, ac$binary_status),
  'prior sact' = rep(NA, 2),
  table(ac$sact_over6weeks_before_studyentry, ac$binary_status),
  'sact during study' = rep(NA, 2),
  table(ac$sact_during_study, ac$binary_status),
  'surgery within 6weeks study' = rep(NA, 2),
  table(ac$surgery_during_study_period, ac$binary_status)
)

# convert to a dataframe ; remove the row names because it corrupts the naming with special characters
summary_airlock <- as.data.frame(matrix, row.names=F) %>% 
  # use the matrix rownames to populate the characteristic column and move this to the start of the table
  mutate(characteristic = rownames(matrix)) %>% relocate(characteristic) %>%
  # rename the columns
  setnames(., colnames(.), c('characteristic', 'no_vte', 'vte')) %>%
  # calculate totals and column proportions
  mutate(overall = no_vte + vte,
         percent_no_vte= round(no_vte/n_no_vte, digits=3),
         percent_vte = round(vte/n_vte, digits=3),
         percent_overall = round(overall/total_n, digits=3)) %>%
  # convert proportions to percentages
  mutate_at(vars(starts_with('percent')), ~round(.*100, 1)) %>%
  # mask any values < 5 as per airlock protocol
  mutate_at(vars(vte, no_vte, overall),
            ~ if_else(. < 5, "<5", as.character(.))
  ) %>%
  # NB need to also mask the counts of any other values which could be used to derive the masked n<5 counts
  # code for this has been redacted as it contained hard-coded values
  # add percentage sign to the percent columns (except for the ones with NA vals)
  mutate_at(vars(starts_with('percent')), ~ if_else(!is.na(.), paste0(., '%'), NA))

# save.
summary_airlock %>% write.table(., "./results/cohortsummary_bystatus.tsv", row.names = FALSE, quote = F, sep = '\t')

# continuous variables
cont_vars <- function(x) { x %>%
  # calculate median, range and IQR in each status group
  summarise(
    # age
    median_age_study_entry = median(age_study_entry),
    range_age_study_entry = paste0(range(age_study_entry), collapse = " - "),
    IQR_age_study_entry = paste0(quantile(age_study_entry, 0.25), '-', quantile(age_study_entry, 0.75)),
    # delay from diagnosis to biopsy
    median_time_diag_to_biopsy = median(time_diag_to_biopsy),
    range_time_diag_to_biopsy = paste0(range(time_diag_to_biopsy), collapse = " - "),
    IQR_time_diag_to_biopsy = paste0(quantile(time_diag_to_biopsy, 0.25), '-', quantile(time_diag_to_biopsy, 0.75)),
    # time to sact
    median_time_to_current_sact = median(time_to_current_sact, na.rm = T),
    range_time_to_current_sact = paste0(range(time_to_current_sact, na.rm = T), collapse = " - "),
    IQR_time_to_current_sact = paste0(quantile(time_to_current_sact, 0.25, na.rm = T), '-', quantile(time_to_current_sact, 0.75, na.rm = T)),
    # time to surgery
    median_time_to_surgery = median(time_to_surgery, na.rm = T),
    range_time_to_surgery = paste0(range(time_to_surgery, na.rm = T), collapse = " - "),
    IQR_time_to_surgery = paste0(quantile(time_to_surgery, 0.25, na.rm = T), '-', quantile(time_to_surgery, 0.75, na.rm = T)),
    # time under observation
    median_time_to_status_months = median(time_to_status_months, na.rm = T),
    range_time_to_status_months = paste0(range(time_to_status_months, na.rm = T), collapse = " - "),
    IQR_time_to_status_months= paste0(quantile(time_to_status_months, 0.25, na.rm = T), '-', quantile(time_to_status_months, 0.75, na.rm = T))
  ) %>%
  # range is a character variable so convert all variables to character to faciliate the pivot
  mutate(across(starts_with('median'), ~as.character(.))) }

summary_cont_vars_grouped <- ac %>%
  #stratify by status
  group_by(binary_status) %>% cont_vars()

# also do an ungrouped one where they aren't separated on status
summary_cont_vars_ungrouped <- ac %>% cont_vars()
  
# join the grouped and ungrouped estimates into one table
summary_cont_vars2 <- rbind(summary_cont_vars_grouped,
                            summary_cont_vars_ungrouped %>% mutate(binary_status = 'overall')) %>%
  # pivot the table so that it has only two columns
  pivot_longer(-binary_status, names_to = "parameter") %>%
  pivot_wider(names_from = binary_status, values_from = value, names_prefix = "status_")

write.table( summary_cont_vars2, "./results/cohortsummary__continuous_vars.tsv", row.names = FALSE, quote = F, sep = '\t')


# VTE rate over time ######
# Supplementary Figure 1

# plot of VTE HAZARD over time (i.e VTE risk)

# function to calculate number at risk at any given timepoint
RiskSetCount <- function(timeindex, survivaltime) {
  atrisk <- NULL
  for (t in timeindex)
    atrisk <- c(atrisk, sum(survivaltime >= t))
  return(atrisk)}

# Supplementary Figure 1
# This estimation uses the muhaz package which uses kernal smoothing and less sensitive to random noise
png('./results/hazard_function_vte_5yr.png', width=1000, height=700)
# define plot parameters. mar = margin = c(bottom, left, top, right) so that text can be written beneath
par(mfrow=c(1,1),mar = c(7,6,2,2)) # define plot parameters
plot(muhaz(ac$time_to_status_months, ac$binary_status, max.time = 60), xlab = 'time months', ylab = 'h(t) - rate of VTE', main = 'Hazard function: instantaneous rate of VTE over time', xaxt = 'n')
axis(side = 1, at=seq(0,60,6), labels = seq(0,60,6))
# to add the risk set grid
grid <- seq(0,60,by= 6)
mtext('Number in risk set', side =1, line = 4, at=0, font =2)
mtext(RiskSetCount(grid,temp$time_to_status_months), side=1, line=5, at=grid, font=2)
dev.off()

# Tumour subtypes #####
# Supplementary Table 2

disease_subtypes <- data.frame(
  # initialise empty df
  disease_type = NA, disease_subtype = NA
)
for (i in (ac %>% pull(disease_type) %>% unique)) {
  x <- data.frame(disease_type = i,
                  disease_subtype = ac %>% filter(disease_type == i) %>%
                    pull(disease_sub_type) %>% unique)
  disease_subtypes <- rbind(disease_subtypes, x)
}
disease_subtypes <- na.omit(disease_subtypes)

disease_subtypes <- ac %>% group_by(disease_type, disease_sub_type) %>% summarise(                                                         freq = n(), .groups = 'drop') %>%
  mutate(freq = case_when(freq < 5 ~ '<5', .default = as.character(freq)))

icd_codes <- list()
for (i in (ac %>% pull(disease_type) %>% unique)) {
  icd_codes[[i]] <- ac %>% filter(disease_type == i) %>% select(site_coded_3char) %>% table
}

# note for HEAD_NECK tumours and 'OTHER', the subtypes shown in supplementary table 2 actually reflect the numbers in the GEL disease_type column (not disease_sub_type columns); code not shown
fwrite(disease_subtypes, './results/summary_tumour_types_airlock.csv', row.names = F)
