# JTH_manuscript
# script for producing tables and figures for JTH manuscript: 
# associations of tumour somatic mutations with VTE
# ncornish 


# config #####
rm(list=ls())
library(kableExtra)
library(gt)
library(metafor)
library(patchwork)

# Figure 1: STROBE flow chart ####
library(DiagrammeR)
library(svglite)
library(rsvg)
library(DiagrammeRsvg)

# load the values for the diagram (see MyWorkflow.Rmd exclusions section for how these are obtained)
strobe_values1 <- fread('STROBE_diagram_vals.tsv')

# get counts
counts <- strobe_values1 %>%
  select(group, count) %>%
  deframe()

strobe_values <- strobe_values1 %>%
  # to make figure clearer group those with missing genetic PCs ('missing covars') and fail_qc
  rbind(data.frame(
    group = 'fail_qc_plus_missing_covars',
    count = counts['fail_qc'] + counts['missing_or_discrepant_covariates']
  ), fill = T) %>%
  # derive the number which passed genetic QC after the various exclusions
  rbind(data.frame(
    group = 'passed_genetic_qc',
    count = counts['all_consenting'] - counts['fail_qc'] - counts['missing_or_discrepant_covariates'] - counts['tumour_sample_prior_2015']
  ), fill = T) %>%
  # derive the number which passed genetic and clinical QC after the various exclusions
  rbind(data.frame(
    group = 'passed_genetic_and_clinical_qc',
    count = counts['all_consenting'] - 
      counts['fail_qc'] - counts['missing_or_discrepant_covariates'] - 
      counts['tumour_sample_prior_2015'] -
      counts['missing_NCRAS'] - counts['incongruent_diagnosis'] -
      counts['benign']
  ), fill = T) 

## HTML version with left aligned text ####
strobe_flow <- grViz("digraph flowchart {
      splines = ortho
      node[fontname = Helvetica, fontsize = 16, shape = rectangle, 
      fixedsize = T]        
      
      include1 [label=
      <<table border='0' cellborder='0' align = 'left'>
      <tr><td align='left'>Cancer programme participants enrolled 2015 - 2019</td></tr>
      <tr><td align='left'>N = @@1 participants with ongoing consent data release v.19 (Oct 2024)</td></tr>
      </table>>]
      
      exclude2 [shape=none, label=
      <<table border='0' cellborder='0' align = 'left'>
      <tr><td align='left'>Genetic data exclusions:</td></tr>
      <tr><td align='left'>@@2 germline or tumour WGS did not pass quality control</td></tr>
      <tr><td align='left'>@@3 tumour sample collected prior to recruitment period</td></tr>
      <tr><td align='left'>@@4 no genetic principal components calculated</td></tr>
      </table>>]
      
      include3 [label=
      <<table border='0' cellborder='0' align = 'left'>
      <tr><td align='left'>@@5 participants with genetic data meeting quality control thresholds</td></tr>
      </table>>] 
      
      exclude4 [shape=none, label=
      <<table border='0' cellborder='0' align = 'left'>
      <tr><td align='left'>Clinical data exclusions:</td></tr>
      <tr><td align='left'>@@6 missing linked national cancer registry record</td></tr>
      <tr><td align='left'>@@7 discrepant diagnostic information between data sources</td></tr>
      <tr><td align='left'>@@8 benign or uncertain histology</td></tr>
      </table>>]
      
      include4 [label=
      <<table border='0' cellborder='0' align = 'left'>
      <tr><td align='left'>@@9 participants with adequate clinical and genetic data </td></tr>
      </table>>] 
      
      exclude5 [shape=none, label=
      <<table border='0' cellborder='0' align = 'left'>
      <tr><td align='left'>Study entry defined as point of tumour biopsy for WGS</td></tr>
      <tr><td align='left'>@@10 excluded: VTE prior to study entry</td></tr>
      </table>>]
      
      include5 [label=
      <<table border='0' cellborder='0' align = 'left'>
      <tr><td align='left'>@@11 pan-cancer participants with no prior history of VTE</td></tr>
      </table>>] 

# create invisible nodes for positioning
      node [shape=none, width=0, height=0, label='']
      p1; p2; p3;
      
      # edge definitions 
      edge [arrowhead=normal, dir=forward, minlen=1]
      
      include1 -> p1 [tailport=sw, headport=nw, arrowhead = none]
      p1 -> exclude2
      p1 -> include3 [tailport=sw, headport=nw]
      
      include3 -> p2 [tailport=sw, headport=nw, arrowhead = none]
      p2 -> exclude4
      p2 -> include4 [tailport=sw, headport=nw]
      
      include4 -> p3 [tailport=sw, headport=nw, arrowhead = none]
      p3 -> exclude5
      p3 -> include5 [tailport=sw, headport=nw]
      
      # force alignment - exclude boxes indented
      {rank=same; include1}
      {rank=same; p1; exclude2} 
      {rank=same; include3}
      {rank=same; p2; exclude4}
      {rank=same; include4}
      {rank=same; p3; exclude5}
      {rank=same; include5}
      }

  # these are R environment variables which must be vectors     
  [1]: strobe_values %>% filter(group == 'all_consenting') %>% dplyr::select(count)
  [2]: strobe_values %>% filter(group == 'fail_qc') %>% dplyr::select(`n (%)`)
  [3]: strobe_values %>% filter(group == 'tumour_sample_prior_2015') %>% dplyr::select(`n (%)`)
  [4]: strobe_values %>% filter(group == 'missing_or_discrepant_covariates') %>% dplyr::select(`n (%)`)
  [5]: strobe_values %>% filter(group == 'passed_genetic_qc') %>% dplyr::select(count)
  [6]: strobe_values %>% filter(group == 'missing_NCRAS') %>% dplyr::select(`n (%)`)
  [7]: strobe_values %>% filter(group == 'incongruent_diagnosis') %>% dplyr::select(`n (%)`)
  [8]: strobe_values %>% filter(group == 'benign') %>% dplyr::select(`n (%)`)
  [9]: strobe_values %>% filter(group == 'passed_genetic_and_clinical_qc') %>% dplyr::select(count)
  [10]: strobe_values %>% filter(group == 'prior_vte') %>% dplyr::select(`n (%)`)
  [11]: strobe_values %>% filter(group == 'primary_cohort') %>% dplyr::select(count)
       ")

# to view image
print(strobe_flow)
# to save
export_svg(strobe_flow) %>% charToRaw %>% rsvg %>% png::writePNG('cancer_cohort_strobe.png')
rm(list=ls())

# Table1: Cohort chracteristics stratified by VTE status at end of follow up ####

# categorical variables
cohort_summary <- fread('cohortsummary_bystatus.tsv') %>% as.data.frame()
cohort_summary_cnt <- fread('cohortsummary_continuous_vars.tsv') %>% as.data.frame()

# SuppTable 4: Univariable Cox PH models ####
clinical <- fread('cox_clinical_univar.tsv') 

clinical2 <- clinical %>%
  # filter for relevant rows
  filter(covariate != 'age_study_entry_scaled' & covariate != 'age_cancerdiag' & covariate != '1'
         & covariate != '2' & covariate != '3' & covariate != '4' & covariate != 'unknown stage') %>%
  # group by term
  group_by(category) %>%
  # insert blank spaces in between category
  reframe(across(everything(), ~ c(NA, .))) %>% ungroup() %>%
  # insert the category name in the blank rows
  mutate(covariate = case_when(is.na(covariate) ~ category, .default= covariate)) %>%
  # change to lower case
  mutate(covariate = tolower(covariate)) %>%
  # remove unecessary characters
  mutate(covariate = gsub('_', ' ', covariate)) %>%
  # remove unecessary characters
  mutate(covariate = gsub('_', ' ', covariate)) %>%
  # remove unecessary characters
  mutate(covariate = gsub('\\(\\)', '', covariate)) %>%
  # round numbers 
  mutate(across(c(HR, l95ci, u95ci), ~format(round(., digits = 2), nsmall = 2))) %>%
  mutate(pval = formatC(pval, format = 'e', digits =2)) %>% 
  mutate('HR[95% CI]' = paste0(HR, '[', l95ci, '-', u95ci, ']')) %>%
  # remove unneccessary cilumns
  dplyr::select(-c(category, model_covars)) %>%
  # rearrange table so sex appears below age
  dplyr::slice(1:2, 25:26, 3:24, 27:nrow(.)) %>%
  # convert NA to blanks
  mutate(across(everything(), ~case_when(str_detect(., 'NA') ~ '', .default = .)))
  
write.csv(clinical2, 'table2_cox_clinical_univar.csv', row.names = F, quote = T)


# Supp Table 5: Gene results #####
## read gene results 
genes_minadj <- fread('cox_genes_minadj.tsv')
genes_fulladj <- fread('cox_genes_fulladj.tsv')
freqs <- fread('frequency_combined_top_genes.tsv')

# first combine minadj and fulladj results
genes_all <- left_join(genes_minadj, genes_fulladj, 
                       by = 'gene', suffix = c('_minadj', '_fulladj')) %>%
  # arrange in order of pval in minadj analysis
  arrange(pval_minadj) %>%
    # format pval as scientific notation
  mutate(across(starts_with('pval'), 
                ~formatC(., format='e', digits=2))) %>%
  # round fdr_p to 4 decimals
  # if fdr_p < 0.0001 then write as this (convert class character)
  mutate(fdr_p = ifelse(fdr_p < 0.0001, '<0.0001', format(round(fdr_p, digits=4), nsmall = 2)),
      # round HR and CIs to decimals
         across(c(starts_with('HR'),starts_with('l95ci'),starts_with('u95ci')),
                ~format(round(., digits=2), nsmall =2))) %>%
  # put model covar columns at the end
  relocate(model_covars_minadj, .before = HR_fulladj) %>% 
  # in model covars put the actual name of each gene
  rowwise() %>%
  mutate(across(starts_with('model_covars'), ~gsub('gene', gene, .))) %>%
  ungroup %>%
  # join with annotations
  left_join(., gene_list_annotated, by = 'gene') %>%
  # clean up the anontation column
  mutate(previously_reported_VTE_locus = case_when(annotation == 'previously reported in VTE GWAS' ~ 
                                                     'previously reported germline VTE locus in GWAS studies',
                                                   annotation == 'previously reported somatic CAT locus' ~
                                                     'previosly reported somatic VTE locus in cancer studies',
                                                   .default = NA)) %>%
  select(-c(annotation)) %>%
  # add a couple of extra columns to facilitate excel highlighting of statistically significant results
  mutate(fdr_p_below_0.1 = case_when(fdr_p < 0.1 ~ T, .default = F),
         p_less_0.05 = case_when(as.numeric(pval_minadj) < 0.05 & 
                                   as.numeric(pval_fulladj) < 0.05 ~ T,
                                              .default = F)) %>%
  relocate(gene_ID, .after = gene)

write.csv(genes_all, 'genes_all.csv', row.names = F)

genes_all <- fread('genes_all.csv')

# Fig 2: forest plots of genes with fdr<0.1#####
genes_minadj <- fread('cox_genes_minadj.tsv')
genes_fulladj <- fread('cox_genes_fulladj.tsv')
freqs <- fread('frequency_combined_top_genes.tsv')
# plot only genes with FDR-P < 0.1
genes_minadj_fdr_sig <- genes_minadj %>% 
  mutate(model = 'min. adj') %>% filter(fdr_p < 0.1) %>%
  # join with mut freqs
  left_join(., freqs %>% filter(cancer == 'PANCANCER'), by = 'gene')
  

gene_factor <- freqs %>% filter(cancer == 'PANCANCER' & 
                                  gene %in% genes_minadj_fdr_sig$gene) %>%
  arrange(desc(prop_mutated)) %>% pull(gene) %>% factor(., levels = unique(.))

df_plot <- rbind(genes_minadj_fdr_sig, 
                # combined minadj and fulladj resutls into single df 
                 genes_fulladj %>% mutate(model = 'full. adj') %>% 
                  filter(gene %in% genes_minadj_fdr_sig$gene) 
                # fill empty columns with NA (FDR-P and mut freq)
                  , fill = T) %>%
  # arrange in order of gene factor
  mutate(gene = factor(gene, levels = gene_factor)) %>% 
  arrange(gene) %>%
  # format pval and prop mutated columns
  mutate(pval = formatC(pval, format = 'e', digits = 2)) %>%
  mutate(prop_mutated = paste0(round(prop_mutated*100), '%')) %>%
  mutate(across(where(is.character), ~case_when(
    is.na(as.character(.)) ~ '',
    . == 'NA%' ~ '',
    .default = .
  ))) %>%
  # group by gene
  group_by(gene) %>%
  # insert blank spaces in between each gene name for ease of viewing on the plot
  reframe(across(everything(), ~ c(., NA))) %>% ungroup() %>%
  # remove the last NA row
  dplyr::slice(1:n()-1)

# set the na action so that it plots the blank rows between categories rather than omitting them
colors <- rep(c("black", "#1F4E79", 
                # pink is a filler colour for the blank rows; i.e. the gaps between the genes
                # will not be used
                "pink"), length(gene_factor)) %>% 
  # delete the last black value
  .[-length(.)]

jpeg('genes_pancancer_forest.jpg', width = 1100, height = 800)
options(na.action = 'na.pass')
forest(x=df_plot$HR, ci.lb = df_plot$l95ci, ci.ub = df_plot$u95ci, 
       slab = df_plot$model, # label on left side
       cex = 2.0, # text size
       font.lab =2,
       refline = 1.0, # line of null effect
       lwd = 2,
       header= c('', 'HR [95% CI]'),
       textpos=c(0, 3.2), # position of annotations
       # extra annotations
       ilab = cbind(df_plot$pval, # pval
                    df_plot$prop_mutated), # mut_freq (%)
       ilab.xpos = c(3.7, -0.3), # position for extra annotations
       pch = 15, # shape of point
       psize = 2.0, # fixed point size
       xlab = 'Hazard ratio for VTE',
       xlim = c(-1.5, 4.0),
       col=colors)

# to add label for P value column
# the nrow(df_plot)+1 gives the y coordinates; font = 2 is bold, cex=1.0 is font size
text(c(3.7, -0.3), nrow(df_plot)+2, c('P', 'Mutation \n Freq (%)'), cex=2.0, font = 2) 
# create category labels
category_labs <- df_plot %>% pull(gene) %>% as.character %>% unique %>% 
  # add the heading 'variable' to the category labels
  c("Gene",.)
# plot category labels on graph
# I used the code below to dynamically figure out the positions for the labels 
cat_positions <- c((nrow(df_plot)+2), # this is the position of the title
                   # then the positions for the category labels
                   seq(from = nrow(df_plot), by = -3, length.out = (length(category_labs)-1)))
for (i in seq_along(cat_positions)) {
  # 0 is the x axis position and pos = 4 aligns the text to the right of the specified position
  # font = 2 is bold
  text(-1.3, cat_positions[i], labels = category_labs[i], pos = 4, font = 2, cex = 2.0)
}
dev.off()
# reset the na action back to default
options(na.action = 'na.omit')


# Supp Figure 3 cuminc plots#########
tp53 <- fread('cuminc_plot_TP53.tsv')
kras <- fread('cuminc_plot_KRAS.tsv')
pcdh15 <- fread('cuminc_plot_PCDH15.tsv')
cdkn2a <- fread('cuminc_plot_CDKN2A.tsv')

fg_res <- fread('FG_genes_minadj.tsv')
fg_res_full <- fread('FG_genes_fulladj.tsv')

# create cuminc function
plot_cum_inc <- function(surv_df, gene_name) {
  # get the SHR and CIs from the fine gray results
  SHR <- fg_res %>% filter(gene == gene_name) %>% mutate(SHR = format(round(SHR, digits =2), nsmall =2)) %>% pull(SHR)
  pval <- fg_res %>% filter(gene == gene_name) %>% mutate(pval = formatC(pval, format = 'e', digits = 2)) %>% pull(pval)
  # combine into an annotation to put on the graph
  annotation <- paste0('SHR = ', SHR, ',\n', 'p = ', pval)
  # now create the plot
  surv_df %>%
    mutate(cum_inc = surv, 
           # I already inverted surv (to make it cuminc) whne creating surv_df
           # but I forgot to do the same for the CIs therefore do it here
           cum_inc_upper = 1 - lower,
           cum_inc_lower = 1 - upper) %>%
    ggplot(aes(x = time, y = cum_inc, color = strata, fill = strata)) +
    scale_x_continuous(breaks = seq(0, 60, 12)) +
    geom_vline(xintercept = 6, color = "gray60", linetype = "dotted", linewidth = 0.8) +
    geom_step() +
    geom_ribbon(aes(ymin = cum_inc_lower, ymax = cum_inc_upper), 
                alpha = 0.3, color = NA) +
    labs(title = gene_name,
         x = "Time (months)",
         y = "VTE cumulative incidence") +
    ylim(0, 0.17) +
    theme_minimal() +
    theme(text = element_text(size = 16),
    axis.text = element_text(size = 16)) +
    guides(color = guide_legend(""), fill = guide_legend("")) #+
}

p1 <- plot_cum_inc(tp53, 'TP53')
p2 <-plot_cum_inc(kras, 'KRAS')
p3 <-plot_cum_inc(cdkn2a, 'CDKN2A')
p4 <-plot_cum_inc(pcdh15, 'PCDH15')

composite <- wrap_plots(p1, p2, p3, p4)
png('genes_cuminc_composit.png', width = 1000, height =800)
composite
dev.off()

# Figure 3 Tumour stratified analyses for each gene ######
cancer_factor <- fread("cohortsummary_bystatus.tsv") %>%
  # pull just the rows relating to cancer
  dplyr::slice(17:32) %>% 
  # remove the category 'OTHER' since this wasn't included in the tumour-stratified analysis
  filter(characteristic != 'OTHER') %>%
  # arrange in order of column 2 (frequencies)
  arrange(desc(`overall, n = N`)) %>% 
  # change formatting of cancer column to remove hyphen
  mutate(characteristic = gsub('_', ' ', characteristic),
         characteristic = gsub('HEAD NECK', 'HEAD AND NECK', characteristic),
         characteristic = gsub('ENDOMETRIAL_CARCINOMA', 'ENDOMETRIAL', characteristic)) %>%
  # convert to a factor keeping the order they appear in this table as
  pull(characteristic) %>% factor(., levels = unique(.))

sample_sizes<- fread('cohortsummary_bystatus.tsv') %>% 
  as.data.frame() %>% slice(17:31) %>% select(c(1,4)) 
colnames(sample_sizes) <- c('cancer', 'n')
# N = total cohort size
pancan <- data.frame(cancer = 'PANCANCER', n = N)
sample_sizes2 <- rbind(sample_sizes, pancan)

freqs <- fread('frequency_combined_top_genes.tsv')

freqs2 <- left_join(freqs, sample_sizes2, by = 'cancer') %>%
  # calculate number of samples with each signature in each caner group
  mutate(n_mutated = round(prop_mutated*n)) %>%
  # remove other as they are not needed for the factoring
  filter(cancer != 'OTHER') %>% select(cancer, gene, n_mutated, prop_mutated)
rm(freqs, sample_sizes, pancan, sample_sizes2)

# read in tumour stratified results
stratified <- fread('cox_top_genes_tumour_stratified.tsv')
df_plot <- stratified %>%
  left_join(., freqs2, by = c('gene', 'cancer')) %>%
  # don't plot any cancers with less than 5 n_mutated
  mutate(across(c(HR, l95ci, u95ci, pval), ~case_when(n_mutated < 5 ~ NA, .default = .))) %>%
  # remove the underscore from cancer names
mutate(cancer = gsub('_', ' ', cancer),
       cancer = gsub('HEAD NECK', 'HEAD AND NECK', cancer),
       cancer = gsub('ENDOMETRIAL_CARCINOMA', 'ENDOMETRIAL', cancer)) %>%
  # round lower CI to 2 digits
  mutate(l95ci = round(l95ci, digits = 2)) %>%
  mutate(gene = factor(gene, levels = c('TP53', 'KRAS', 'PCDH15', 'CDKN2A' )),
         cancer = factor(cancer, levels = cancer_factor)) %>%
  # Arrange by cancer and gene for consistent ordering
  arrange(cancer, gene) %>%
  mutate(
# set positions for each cancer  
    cancer_pos = (length(levels(cancer)) - as.numeric(cancer)) * 2 + 1,
    gene_offset = (length(levels(gene)) - as.numeric(gene) + 1 - mean(length(levels(gene)) - as.numeric(gene) + 1)) * 0.3,
    # Add alternating band indicator
    band_fill = ifelse(as.numeric(cancer) %% 2 == 0, "gray95", "white")
  ) 

p <- df_plot %>% 
           ggplot(aes(x = HR, y = cancer_pos + gene_offset, color = gene)) +
  # Add alternating background using geom_tile
  geom_tile(aes(x=1, y = cancer_pos, fill = band_fill
                ), 
            width = Inf, height = 1.5, alpha = 0.3) +
  scale_fill_identity() +
  # Add vertical line at HR = 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  # Add confidence intervals
  geom_errorbarh(aes(xmin = l95ci, xmax = u95ci), height = 0.1, alpha = 0.8) +
  # Add HR points
  geom_point(size = 3, alpha = 0.9) +
  # Set x-axis to log2 scale (doubling scale)
  scale_x_continuous(
    trans = "log2",
    breaks = c(0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16),
    labels = c("0.0625", "0.125", "0.25", "0.5", "1", "2", "4", "8", "16")
  ) +
  # Set y-axis with cancer names and larger spacing
  scale_y_continuous(
    breaks = seq(1, length(levels(df_plot$cancer)) * 2 - 1, by = 2),
    labels = rev(levels(df_plot$cancer))
  ) +
  # Color scale for genes
  scale_color_manual(values = c("KRAS" = "#1F4E79", 
                              "CDKN2A" = "#000",
                                "TP53" = "#9E2A2B",
                                "PCDH15" = "#4E6E58")) +
  # Labels and theme
  labs(
    x = "Hazard Ratio for VTE",
    y = "",
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18, face = "bold"),
    legend.title= element_text(size = 18, face = "bold"),
    legend.position = "bottom",
    legend.key = element_rect(fill = NA, color = NA))

ggsave("genes_tumour_strat.jpg", plot = p, width = 12, height = 12, dpi = 300)

# Supp Table 6 Gene sensitivity analyses ##########
# full adj pan cancer
df1 <- fread('cox_genes_minadj.tsv') %>%
  mutate(sensitivity_analysis = 'min.adj pan-cancer') %>%
  filter(gene %in% c('TP53', 'KRAS', 'CDKN2A', 'PCDH15')) %>%
  mutate(stratum = 'whole cohort')
df2 <- fread('cox_genes_fulladj.tsv') %>%
  mutate(sensitivity_analysis = 'full.adj pan-cancer') %>%
  filter(gene %in% c('TP53', 'KRAS', 'CDKN2A', 'PCDH15')) %>%
  mutate(stratum = 'whole cohort')
# pancancer time stratified
df3 <- fread('cox_genes_time_stratified.tsv') %>%
  mutate(sensitivity_analysis = 'time-stratified') %>% dplyr::rename('gene' =term) %>%
  mutate(stratum = paste('year', time_point))
# pancancer new untreated cancer (<42 days diagnosis to study entry)
# this file has now been updated
df4_all <- fread('cox_new_untreated_timezero_diag_genes.csv') %>%
  mutate(sensitivity_analysis = 'min.adj time from cancer diagnosis to VTE',
         stratum = 'delay < 42 days from diagnosis to study and no prior treatment') %>%
  select(-fdr_p)
df4_all_adj <- fread('cox_new_untreated_timezero_diag_genes_fulladj.csv') %>%
        mutate(stratum = 'delay < 42 days from diagnosis to study and no prior treatment',
               sensitivity_analysis = 'full.adj from cancer diagnosis to VTE'
                )

df4 <- rbind(df4_all, df4_all_adj) %>% filter(gene %in% c('TP53', 'KRAS', 'CDKN2A', 'PCDH15'))
rm(df4_all, df4_all_adj)
# no indication for anticoag
df5 <- fread('cox_genes_anticoag_stratified.tsv') %>%
  mutate(sensitivity_analysis = 'possible history of anticoagulant or antiplatelet therapy') %>%
  mutate(stratum = case_when(prior_anticoag_indication == 'yes' ~ 'medical indication for anticoagulation or antiplatelets',
                             prior_anticoag_indication == 'no' ~ 'no medical indication for anticoagulation or antiplatelets'))

# ancestry stratified
df6 <- fread('cox_genes_ancestry_stratified.tsv') %>%
  mutate(sensitivity_analysis = 'ancestry-stratified') %>%
  mutate(stratum = case_when(ancestry == 'European' ~ 'European ancestry only',
                             ancestry == 'African' ~ 'African ancestry only',
                             ancestry == 'South Asian' ~ 'South Asian ancestry only',
                             ancestry == 'East Asian' ~ 'East Asian ancestry only',
                             ancestry == 'American' ~ 'American ancestry only')) %>%
  # for East Asians regression models did not converge for any gene due to data sparsity so remove from the table,
  filter(stratum != 'East Asian ancestry only')

# multi-gene regression
df7 <- fread('cox_multigene_regression.tsv') %>%
  mutate(sensitivity_analysis = 'multi-gene regression') %>%
  mutate(stratum = 'whole cohort')
# tumour stratified regression
df8 <- fread('cox_top_genes_tumour_stratified.tsv') %>%
  mutate(sensitivity_analysis = 'tumour-stratified') %>%
  mutate(stratum = cancer)
# fine and gray
df9 <- fread('FG_genes_minadj.tsv') %>%
  mutate(sensitivity_analysis = 'Fine and Gray minimally adjusted model') %>%
  mutate(stratum = 'whole cohort') %>%
  dplyr::rename('model_covars'= model)
df10 <- fread('FG_genes_fulladj.tsv') %>%
  mutate(sensitivity_analysis = 'Fine and Gray fully adjusted model') %>%
  mutate(stratum = 'whole cohort') %>%
  dplyr::rename('model_covars'= model)

df <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, fill = T) %>%
  # combine all tables, concatenate HR and CI columns and round to 2 digits
  mutate(across(c(HR, l95ci, u95ci), ~format(round(., digits=2), nsmall =2))) %>%
  mutate(pval = formatC(pval, format='e', digits =2)) %>%
  mutate('HR [95% CI]' = paste0(HR, ' [', l95ci, ',', u95ci, ']')) %>%
  # select relevant columns
  dplyr::select(gene, sensitivity_analysis, stratum, model_covars, 'HR [95% CI]', pval) %>%
  # standardise the model covars column so it has 'gene' instead of the gene name 
  # otherwise the column has to be duplicated for each gene
  mutate(model_covars = gsub('TP53|KRAS|CDKN2A|PCDH15', 'gene', model_covars)) %>%
  mutate(model_covars = gsub('pc1 \\+ pc2 \\+ pc3 \\+ pc4', 'pc1-4', model_covars)) %>%
  # pivot wider so each gene appears on only one row
  pivot_wider(
    names_from = gene,
    values_from = c('HR [95% CI]', pval),
    names_sep = "_"
  ) %>%
  relocate(pval_TP53, .after = `HR [95% CI]_TP53`) %>%
  relocate(pval_PCDH15, .after = `HR [95% CI]_PCDH15`) %>%
  relocate(pval_CDKN2A, .after = `HR [95% CI]_CDKN2A`) %>%
  relocate(pval_KRAS, .after = `HR [95% CI]_KRAS`) %>%
  # replace NA values with blank for cleaner output
  mutate(across(c(4:ncol(.)), ~case_when(str_detect(., 'NA') ~ '', .default = .))) %>%
  relocate(model_covars, .after = last_col())

fwrite(df, 'gene_sensitivity_analyses.csv', row.names = F)
rm(df1, df2, df3, df4, df5, df6, df7, df8, df9)

# ANOVA and interaction P-values from tumour*gene interaction tests 
anovas <- fread('anova_tumour_gene_interactions.tsv') 
inter <- fread('interaction_terms_tumour_gene.tsv') 

# Supp Table 7 TMB main and sensitivity analyses ##########
# min adj pan cancer
df1 <- fread('cox_TMB_minadj.tsv') %>%
  mutate(sensitivity_analysis = 'min.adj pan-cancer') %>%
  filter(str_detect(term, 'tmb_categories')|str_detect(term, 'tmb_binary_thresh20tmb_above_20')) %>%
  mutate(stratum = 'whole cohort')
# full adj pan cancer
df2 <- fread('cox_TMB_fulladj.tsv') %>%
  mutate(sensitivity_analysis = 'full.adj pan-cancer') %>%
  filter(str_detect(term, 'tmb_categories')|str_detect(term, 'tmb_binary_thresh20tmb_above_20')) %>%
  mutate(stratum = 'whole cohort') %>%
  # add a model covars column
  mutate(model_covars = 'tmb + age_study_entry + sex + pc1-4 + tumour_type + stage_grouped + current_sact + sact_over6weeks_before_studyentry')
# pancancer time stratified
df3 <- fread('cox_TMB_time_stratified.tsv') %>%
  mutate(sensitivity_analysis = 'time-stratified') %>% 
  filter(str_detect(term, 'tmb_categories')|str_detect(term, 'tmb_binary_thresh20tmb_above_20')) %>%
  mutate(stratum = paste('year', time_point))

# pancancer new untreated cancer (<42 days diagnosis to study entry)
df4 <- fread('cox_new_untreated_timezero_diag_TMB.tsv') %>%
  mutate(sensitivity_analysis = 'time zero at cancer diagnosis') %>%
  mutate(stratum = 'delay < 42 days from diagnosis to study and no prior treatment')

# no indication for anticoag
df5 <- fread('cox_TMB_anticoag_stratified.tsv') %>%
  filter(str_detect(term, 'tmb_categories')|str_detect(term, 'tmb_binary_thresh20tmb_above_20')) %>%
  mutate(sensitivity_analysis = 'possible history of anticoagulant or antiplatelet therapy') %>%
  mutate(stratum = case_when(prior_anticoag_indication == 'yes' ~ 'medical indication for anticoagulation or antiplatelets',
                             prior_anticoag_indication == 'no' ~ 'no medical indication for anticoagulation or antiplatelets'))

# ancestry stratified
df6 <- fread('cox_TMB_ancestry_stratified.tsv') %>%
  filter(str_detect(term, 'tmb_categories')|str_detect(term, 'tmb_binary_thresh20tmb_above_20')) %>%
  mutate(sensitivity_analysis = 'ancestry-stratified') %>%
  mutate(stratum = case_when(ancestry == 'European' ~ 'European ancestry only',
                             ancestry == 'African' ~ 'African ancestry only',
                             ancestry == 'South Asian' ~ 'South Asian ancestry only',
                             ancestry == 'East Asian' ~ 'East Asian ancestry only',
                             ancestry == 'American' ~ 'American ancestry only')) %>%
  # for East Asians regression models did not convergedue to data sparsity so remove from the table,
  filter(stratum != 'East Asian ancestry only')
# df7 ICI interactions
df7 <- fread('cox_TMB_ICI_interactions.tsv') %>%
  mutate(sensitivity_analysis = 'adjustment for immune check point inhibitors') %>%
  # filter just for the model where I have included ICI
  filter(model == 'tmb_binary + ICI + tmb_binary_thresh20:ICI' &
           term == 'tmb_binary_thresh20tmb_above_20') %>%
  mutate(stratum = 'whole cohort') %>% dplyr::rename('model_covars' = model)
# tumour stratified regression
df8 <- fread('cox_TMB_tumour_stratified.tsv') %>%
  filter(str_detect(term, 'tmb_categories')|str_detect(term, 'tmb_binary_thresh20tmb_above_20')) %>%
  mutate(sensitivity_analysis = 'tumour-stratified') %>%
  mutate(stratum = cancer) %>%
  mutate(across(c(HR, l95ci, u95ci, pval), 
                # for these cancer types, samples with tmb>20 were very rare (n < 5) and regression models did not converge
                ~case_when(cancer == 'HEPATOPANCREATOBILIARY'|
                             cancer == 'RENAL'|
                             cancer == 'SARCOMA'|
                             cancer == 'PROSTATE'|
                             cancer == 'HEAD_NECK'|
                             cancer == 'HAEM_ONC' ~ NA_real_, 
                           .default = .)))
# fine and gray
df9 <- fread('FG_TMB_minadj.tsv') %>%
  mutate(sensitivity_analysis = 'Fine and Gray minimally adjusted model') %>%
  mutate(stratum = 'whole cohort') %>%
  dplyr::rename('model_covars'= model) %>%
  filter(term == 'tmb_binary_thresh20tmb_above_20') %>% 
  # change SHR to HR to facilitate rbind
  rename('HR' = SHR)
df10 <- fread('FG_TMB_fulladj.tsv') %>%
  mutate(sensitivity_analysis = 'Fine and Gray fully adjusted model') %>%
  mutate(stratum = 'whole cohort') %>%
  dplyr::rename('model_covars'= model) %>%
   filter(term == 'tmb_binary_thresh20tmb_above_20') %>% 
  # change SHR to HR to facilitate rbind
  rename('HR' = SHR)

df <- rbind(df1, df2, df7, df3, df4, df5, df6, df8, df9, df10, fill = T) %>%
  # combine all tables, concatenate HR and CI columns and round to 2 digits
  mutate(across(c(HR, l95ci, u95ci), ~format(round(., digits=2), nsmall =2))) %>%
  mutate(pval = formatC(pval, format='e', digits =2)) %>%
  mutate('HR [95% CI]' = paste0(HR, ' [', l95ci, ',', u95ci, ']')) %>%
  # select relevant columns
  dplyr::select(term, sensitivity_analysis, stratum, model_covars, 'HR [95% CI]', pval) %>%
  # standardise the model covars column so it has 'gene' instead of the gene name 
  # otherwise the column has to be duplicated for each gene
  mutate(model_covars = gsub('tmb_binary_thresh20|tmb_categories|tmb_binary', 'TMB', model_covars),
         model_covars = gsub('pc1 \\+ pc2 \\+ pc3 \\+ pc4', 'pc1-4', model_covars),
         # clen the term names
         term = case_when(term == 'tmb_binary_thresh20tmb_above_20' ~ 'TMB>=20_binary',
                          term == 'tmb_categoriestmb_between_5_10' ~ 'TMB 5-9',
                          term == 'tmb_categoriestmb_between_10_20' ~ 'TMB 10-19',
                          term == 'tmb_categoriestmb_above_20' ~ 'TMB >=20')) %>%
  # pivot wider so each term appears on only one row
  pivot_wider(
    names_from = term,
    values_from = c('HR [95% CI]', pval),
    names_sep = "_"
  ) %>% 
  # re-order columns
  relocate(`HR [95% CI]_TMB>=20_binary`, `pval_TMB>=20_binary`,
           `HR [95% CI]_TMB >=20`, `pval_TMB >=20`,
           `HR [95% CI]_TMB 10-19`, `pval_TMB 10-19`,
           `HR [95% CI]_TMB 5-9`, `pval_TMB 5-9`, model_covars, .after = stratum) %>%
  # replace NA values with blank for cleaner output
  mutate(across(c(3:ncol(.)), ~case_when(str_detect(., 'NA') ~ '', .default = .)))


fwrite(df, 'TMB_sensitivity_analyses.csv', row.names = F)

# Supp Figure 4 #####
# Done inside GEL RE with raw data (box plot)

# Supp Figure 5 TMB leave one out analyses##########
tmb_exclude_highTMB_cancers <- read.table('cox_TMB_exclude_highTMB_cancers.tsv')
tmb_summ_stats <- fread('TMB_summarystats.tsv')
tmb_ranks <- tmb_summ_stats %>% arrange(desc(tmb_median)) %>%
  filter(!tumour %in% c('PANCANCER', 'OTHER')) %>%
  pull(tumour) %>% factor(., levels = unique(.))
#Leave one out forest
df_plot <- fread('cox_TMB_leave_one_cancer_out.tsv') %>% 
  filter(!cancer_left_out %in% c('PANCANCER', 'OTHER')) %>%
  # arrange in order of median TMB
  mutate(cancer_left_out =  factor(cancer_left_out, levels = tmb_ranks)) %>%
  arrange(cancer_left_out) 
png('TMB_leave_one_out_forest.png', width = 700, height = 600)
forest(x=df_plot$HR, ci.lb = df_plot$l95ci, ci.ub = df_plot$u95ci, 
       slab = df_plot$cancer_left_out, # label on left side
       textpos=c(-0.3, 1.4), # position of annotations
       cex = 1, # text size
       refline = 1.0, # line of null effect
       #      main = 'Effect of somatic gene mutations on rate of VTE', # title for the overall graph
       header= c('tumour excluded', 'HR [95% CI]'),
       pch = 15, # shape of point
       psize = 2.0, # fixed point size
       xlab = 'Hazard ratio for VTE',
       xlim = c(-0.5, 1.5))
dev.off()

# Supp Table8 signature results ####
sigs_minadj <- fread('cox_signatures_minadj.tsv') %>% mutate(model = 'min.adj')
sigs_fulladj <- fread('cox_signatures_fulladj.tsv')  %>% mutate(model = 'full.adj')

# https://cancer.sanger.ac.uk/signatures/sbs/
# Signature annotations derived from Sanger
sig_annotation <- fread('COSMIC_signatures_aetiology.csv')

df <- rbind(sigs_minadj, sigs_fulladj, fill = T) %>%
         mutate(
         model_covars = gsub('pc1 \\+ pc2 \\+ pc3 \\+ pc4', 'pc1-4', model_covars)) %>%
  # concatenate HR and CI columns and round to 2 digits
  mutate(across(c(HR, l95ci, u95ci), ~format(round(., digits=2), nsmall =2))) %>%
  mutate(pval = formatC(pval, format='e', digits =2)) %>%
  mutate('HR [95% CI]' = paste0(HR, ' [', l95ci, ',', u95ci, ']')) %>%
  select(-c(HR, contains('95ci'))) %>%
    # pivot wider so each signature appears on only one row
  pivot_wider(
    names_from = model,
    values_from = c('HR [95% CI]', pval, model_covars),
    names_sep = "_"
  ) %>%
  # select relevant columns
  select(term, `HR [95% CI]_min.adj`, `pval_min.adj`, model_covars_min.adj,
         `HR [95% CI]_full.adj`, `pval_full.adj`, model_covars_full.adj
         ) %>%
  left_join(., sig_annotation, by = c('term' = 'signature')) %>%
  mutate(term = gsub('signature_', 'SBS', term))
  
signature_order <- c('SBS1','SBS2','SBS3','SBS4','SBS5','SBS6','SBS7',
                     'SBS8','SBS9','SBS10','SBS11','SBS12','SBS13','SBS14',
                     'SBS15','SBS16','SBS17','SBS18','SBS19','SBS20','SBS21','SBS22',
                     'SBS23','SBS24','SBS25','SBS26','SBS27','SBS28','SBS29','SBS30')

df2 <- df %>% mutate(signature = factor(term, levels = signature_order)) %>% 
  mutate(fdr_p_min.adj = round(p.adjust(pval_min.adj, method = 'fdr'), digits =2)) %>%
  select(-c(term)) %>% 
  relocate(signature, signature_aetiology) %>% 
  relocate(fdr_p_min.adj, .after = pval_min.adj) %>%
  arrange(signature) 

fwrite(df2, 'signatures_main_results.csv', row.names = F)

# Supp Table 9 sensitivity analyses ######

# frequency of each signature by tumour type
freqs <- fread('mut_sig_summstats_airlock.tsv')

sample_sizes<- fread('cohortsummary_bystatus.tsv') %>%
  as.data.frame() %>% slice(17:31) %>% select(c(1,4))
colnames(sample_sizes) <- c('cancer', 'n')
# N = total cohort size
pancan <- data.frame(cancer = 'PANCANCER', n = N)
sample_sizes2 <- rbind(sample_sizes, pancan)

freqs2 <- left_join(freqs, sample_sizes2, by = 'cancer') %>%
# calculate number of samples with each signature in each caner group
    mutate(prop_sig = percent_with_sig/100) %>%
  mutate(n_sig = round(prop_sig*n)) %>%
  # remove other as they are not needed for the factoring
  filter(cancer != 'OTHER') %>% select(cancer, signature, n_sig, percent_with_sig) %>% 
  dplyr::rename('term' = signature)
rm(freqs, sample_sizes, pancan, sample_sizes2)
# full adj pan cancer
# FDR sig results
my_sigs <- c('signature_6', 'signature_8', 'signature_19', 'signature_26')
df1 <- fread('cox_signatures_minadj.tsv') %>%
  mutate(sensitivity_analysis = 'min.adj pan-cancer') %>%
  filter(term %in% my_sigs) %>%
  mutate(stratum = 'whole cohort')
df2 <- fread('cox_signatures_fulladj.tsv') %>%
  mutate(sensitivity_analysis = 'full.adj pan-cancer') %>%
  filter(term %in% my_sigs) %>%
  mutate(stratum = 'whole cohort')
# pancancer time stratified
df3 <- fread('cox_signatures_time_stratified.tsv') %>%
  mutate(sensitivity_analysis = 'time-stratified') %>% 
  mutate(stratum = paste('year', time_point))
# pancancer new untreated cancer (<42 days diagnosis to study entry)
df4 <- fread('cox_new_untreated_timezero_diag_signatures.tsv') %>%
  filter(term %in% my_sigs) %>%
  mutate(sensitivity_analysis = 'time zero at cancer diagnosis') %>%
  mutate(stratum = 'delay < 42 days from diagnosis to study and no prior treatment')

# no indication for anticoag
df5 <- fread('cox_signatures_anticoag_stratified.tsv') %>%
  mutate(sensitivity_analysis = 'possible history of anticoagulant or antiplatelet therapy') %>%
  mutate(stratum = case_when(prior_anticoag_indication == 'yes' ~ 'medical indication for anticoagulation or antiplatelets',
                             prior_anticoag_indication == 'no' ~ 'no medical indication for anticoagulation or antiplatelets'))

# ancestry stratified
df6 <- fread('cox_signatures_ancestry_stratified.tsv') %>%
  mutate(sensitivity_analysis = 'ancestry-stratified') %>%
  mutate(stratum = case_when(ancestry == 'European' ~ 'European ancestry only',
                             ancestry == 'African' ~ 'African ancestry only',
                             ancestry == 'South Asian' ~ 'South Asian ancestry only',
                             ancestry == 'East Asian' ~ 'East Asian ancestry only',
                             ancestry == 'American' ~ 'American ancestry only')) %>%
  filter(stratum != 'East Asian ancestry only')
# ICI
df7 <- fread('cox_signature_ICI_interactions.tsv') %>%
  mutate(sensitivity_analysis = 'adjustment for immune check point inhibitors',
         stratum = 'whole cohort') %>%
  # filter for only the model with the time interaction term
  filter(grepl(":ICI", model)) %>% 
  # filter specifically for the signature term
  filter(term %in% my_sigs) %>%
  # clean the model covars column
  mutate(model_covars = gsub("_\\d{1,2}", "", model)) %>% select(-c(model))

# tumour stratified regression
df8 <- fread('cox_signatures_tumour_stratified.tsv') %>%
  mutate(sensitivity_analysis = 'tumour-stratified') %>%
  mutate(stratum = cancer) %>%
  # remove results for any cancers where n_sig < 5 as regression model too unstable to converge
  left_join(., freqs2, by = c('cancer', 'term')) %>% 
  mutate(across(c(HR, u95ci, l95ci, pval), ~case_when(n_sig < 5 ~ NA_real_,
                                                      .default = .))) %>%
  # now remove unnecessary columms
  select(-c(n_sig, percent_with_sig))

df9 <- fread('cox_signature_continuous_nonzeros.tsv') %>%
  mutate(sensitivity_analysis = 'signature contribution to total TMB') %>%
  mutate(stratum = 'signature present')

df10 <- fread('FG_signatures_minadj.tsv') %>%
  mutate(sensitivity_analysis = 'Fine and Gray minimally adjusted model') %>%
  mutate(stratum = 'whole cohort') %>%
  # change SHR to HR to facilitate rbind
  rename('HR' = SHR) %>% 
  mutate(model_covars = 'signature_6 + age_study_entry + sex + pc1-4')
df11 <- fread('FG_signatures_fulladj.tsv') %>%
  mutate(sensitivity_analysis = 'Fine and Gray fully adjusted model') %>%
  mutate(stratum = 'whole cohort') %>%
  # change SHR to HR to facilitate rbind
  rename('HR' = SHR) %>% 
  mutate(model_covars = 'signature_6 + age_study_entry + sex + pc1-4 + disease_type + stage_grouped + current_sact + sact_over6weeks_before_studyentry')

df <- rbind(df1, df2, df7, df9, df3, df4, df5, df6, #df7,
            df8, df10, df11, fill = T) %>%
  # combine all tables, concatenate HR and CI columns and round to 2 digits
  mutate(across(c(HR, l95ci, u95ci), ~format(round(., digits=2), nsmall =2))) %>%
  mutate(pval = formatC(pval, format='e', digits =2)) %>%
  mutate('HR [95% CI]' = paste0(HR, ' [', l95ci, ',', u95ci, ']')) %>%
  # select relevant columns
  dplyr::select(-c(HR, l95ci, u95ci, time_point, ncase, prior_anticoag_indication, ancestry, cancer)) %>%
  relocate(pval, .after = `HR [95% CI]`) %>%
    mutate(
          model_covars = gsub('pc1 \\+ pc2 \\+ pc3 \\+ pc4', 'pc1-4', model_covars),
         signature = gsub('signature_', 'SBS', term)
         ) %>% select(-c(term))

# pivot wider so each term appears on only one row
df_new <- df %>%  pivot_wider(
    names_from = signature,
    values_from = c('HR [95% CI]', pval),
    names_sep = "_"
  ) %>%
  # replace NA values with blank for cleaner output
   select(sensitivity_analysis, stratum,
          `HR [95% CI]_SBS6`, pval_SBS6,
          `HR [95% CI]_SBS8`, pval_SBS8,
          `HR [95% CI]_SBS19`, pval_SBS19,
          `HR [95% CI]_SBS26`, pval_SBS26,
          model_covars)
  
fwrite(df_new, 'signature_sensitivity_analyses.csv', row.names = F)

# have only included time stratified analyses. These results show the time dependent interaction term for TP53
time_int <- fread('cox_signatures_time_interaction.tsv') 
# no significant time interactions
time_int %>% arrange(term) %>%
  filter(pval<0.05)




# Table 2 combined sig and TMB #####
df1 <- fread('cox_TMB_minadj.tsv') %>%
  filter(str_detect(term, 'tmb_categories')|str_detect(term, 'tmb_binary_thresh20tmb_above_20')) %>% 
  # round numbers 
  mutate(across(c(HR, l95ci, u95ci), ~format(round(., digits = 2), nsmall = 2))) %>%
  mutate(pval = format(round(pval, digits=3), nsmall=3)) %>% 
  mutate('HR[95% CI], P' = paste0(HR, '[', l95ci, '-', u95ci, '], P=', pval)) %>%
  select(term, `HR[95% CI], P`)

df2 <- fread('cox_TMB_fulladj.tsv') %>%
  filter(str_detect(term, 'tmb_categories')|str_detect(term, 'tmb_binary_thresh20tmb_above_20')) %>%
  # round numbers 
  mutate(across(c(HR, l95ci, u95ci), ~format(round(., digits = 2), nsmall = 2))) %>%
  mutate(pval = format(round(pval, digits=3), nsmall =3)) %>% 
  mutate('HR[95% CI], P' = paste0(HR, '[', l95ci, '-', u95ci, '], P=', pval)) %>%
  select(term, `HR[95% CI], P`)

df_tmb <- full_join(df1, df2, by = 'term', suffix = c('_min.adj', '_full.adj')) %>% 
  mutate(category = case_when(str_detect(term, 'tmb_categories') ~'tmb_category',
                              str_detect(term, 'tmb_binary_thresh20tmb_above_20') ~ 'tmb_binary'))
rm(df1, df2)

my_sigs <- fread('cox_signatures_minadj.tsv') %>% 
  mutate(fdr_p = p.adjust(pval, method = 'fdr')) %>% filter(fdr_p <0.1) %>%
  pull(term)
sig_annotation <- fread('COSMIC_signatures_aetiology_fromClaude.csv') %>%
  mutate(signature = gsub('signature_', 'SBS', signature))

df1 <- fread('cox_signatures_minadj.tsv') %>%
  filter(term %in% my_sigs) %>% 
  # round numbers 
  mutate(across(c(HR, l95ci, u95ci), ~format(round(., digits = 2), nsmall = 2))) %>%
  mutate(pval = format(round(pval, digits=3), nsmall =3)) %>% 
  mutate('HR[95% CI], P' = paste0(HR, ' [', l95ci, '-', u95ci, '], P=', pval)) %>%
  select(term, `HR[95% CI], P`)
df2 <- fread('cox_signatures_fulladj.tsv') %>%
  filter(term %in% my_sigs) %>% 
  # round numbers 
  mutate(across(c(HR, l95ci, u95ci), ~format(round(., digits = 2), nsmall = 2))) %>%
  mutate(pval = format(round(pval, digits=3), nsmall =3)) %>% 
  mutate('HR[95% CI], P' = paste0(HR, ' [', l95ci, '-', u95ci, '], P=', pval)) %>%
  select(term, `HR[95% CI], P`)
df_sigs <- full_join(df1, df2, by = 'term', suffix = c('_min.adj', '_full.adj')) %>% 
  mutate(category = 'SBS mutational signature') %>%
  mutate(term = factor(gsub('signature_', 'SBS', term), 
                       levels = c('SBS4', 'SBS5', 'SBS6', 
                                  'SBS8', 'SBS17', 'SBS19', 
                                  'SBS26', 'SBS30'))) %>%
  arrange(term) %>% mutate(term = as.character(term)) %>%
  left_join(., sig_annotation, by = c('term'='signature')) %>%
  mutate(term = paste0(term, ' (', signature_aetiology, ')')) %>% select(-signature_aetiology)
rm(df1, df2)           

# merge tables
combined <- rbind(df_tmb, df_sigs) %>%
  group_by(category) %>%
  # insert blank spaces in between category
  reframe(across(everything(), ~ c(NA, .))) %>% ungroup() %>%
  # insert the category name in the blank rows
  mutate(term = case_when(is.na(term) ~ category, .default= term)) %>%
  # convert NA to blanks
  select(-category)
# was originally caleld table 3 but has moved in manuscript
write.csv(combined, 'table3_TMB_sigs.csv', row.names = F, quote = T)

# Supp Figure 6 mutational signatures continuous distribution across tumour types ######
# get the order of the cancers (based on sample size)
cancer_factor <- fread("cohortsummary_bystatus.tsv") %>%
  # pull just the rows relating to cancer
  dplyr::slice(17:32) %>% 
  # remove the category 'OTHER' since this wasn't included in the tumour-stratified analysis
  filter(characteristic != 'OTHER') %>%
  # arrange in order of column 2 (frequencies)
  arrange(desc(`overall, n = N`)) %>% 
  # convert to a factor keeping the order they appear in this table as
  pull(characteristic) %>% c('PANCANCER',.) %>%
  factor(., levels = unique(.))

# Mutsig summstats:
df <- fread('mut_sig_summstats_airlock.tsv') %>%
  mutate(term= gsub('signature_', 'SBS', signature)) %>%
         arrange(term) %>%
  # order cancer as per cancer factor: see below
  mutate(cancer = factor(cancer, levels = cancer_factor),
         percent = percent_with_sig) %>%
filter(!cancer == 'OTHER') %>%
  select(term, cancer, median, min, max, q1, q3, percent)

# TMB summstats. 
df1 <- fread('TMB_summarystats.tsv') %>%
  # order cancer as per cancer factor: see below
  mutate(cancer = factor(tumour, levels = cancer_factor),
         median = tmb_median,
         term = 'TMB',
         # tmb_very_high_prop is a proportion (decimal) not a percentage so need to convert
         percent = tmb_very_high_prop*100,
         min = tmb_min,
         max = tmb_max) %>%
  select(term, cancer, median, min, max, q1, q3, percent)

# load the freqs2 file above  
freqs <- fread('mut_sig_summstats_airlock.tsv')
sample_sizes<- fread('cohortsummary_bystatus.tsv') %>%
  as.data.frame() %>% slice(17:31) %>% select(c(1,4))
colnames(sample_sizes) <- c('cancer', 'n')
# N = total cohort size
pancan <- data.frame(cancer = 'PANCANCER', n = N)
sample_sizes2 <- rbind(sample_sizes, pancan)

freqs2 <- left_join(freqs, sample_sizes2, by = 'cancer') %>%
  # calculate number of samples with each signature in each caner group
  mutate(prop_sig = percent_with_sig/100) %>%
  mutate(n_sig = round(prop_sig*n)) %>%
  # remove other as they are not needed for the factoring
  filter(cancer != 'OTHER') %>% select(cancer, signature, n_sig, percent_with_sig) %>% 
  dplyr::rename('term' = signature)
rm(freqs, sample_sizes, pancan, sample_sizes2)

df_plot <- df %>% left_join(., 
                            freqs2 %>% mutate(term= gsub('signature_', 'SBS', term)) %>% 
                              select(term, cancer, n_sig),
                            by = c('term', 'cancer')) %>% 
  rbind(., df1, fill = T) %>%
  mutate(term = factor(term, levels = c('TMB', 'SBS6', 'SBS8', 'SBS19', 'SBS26'))) %>%
  # if n_sig was < 5 in freqs then make the values NA because I don't want to plot them,
  mutate(across(c(median, q1, q3, min, max, percent), ~case_when(n_sig < 5 ~ NA_real_,
                                                                .default = .)))
rm(df, df1)


# create the plot
p <- df_plot %>% filter(term != 'TMB') %>%
  ggplot(aes(x = 1, y = median)) +
  geom_point(aes(y = median), size = 3, alpha = 0.8) +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 0.3, alpha = 0.8) +
  geom_text(aes(y = q3*1.005, label = paste0(round(percent, 1), "%")), 
            vjust = -0.3, size = 5) +
  facet_grid(term ~ cancer, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14),
    strip.text.x = element_text(size = 14, angle = 90),
    strip.text.y = element_text(size = 20),
    axis.title.y = element_text(size= 20),
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) +
  labs(x = "",
       y = "Contribution of signature to total mutation burden (%)"
       ) 

png('signature_distributions.png', width = 800, height = 1000)
p
dev.off()

# Figure 4 combined TMB and mut sigs stratified by cancer ######
df_plot2 <- df_plot %>% 
  # remove 'OTHER' cancer
filter(cancer != 'OTHER') %>%
mutate(term = case_when(term == 'TMB' ~ 'TMB>=20', .default = term),
                               term = factor(term, levels = c('TMB>=20', 'SBS6', 'SBS19', 'SBS26', 'SBS8')),
                               cancer = factor(cancer, levels = cancer_factor)) %>%
 mutate(percent = case_when(percent == 0 ~ NA_real_, .default = percent)) %>%
  arrange(term, cancer) 
  
p2 <- ggplot(df_plot2, aes(x = cancer, y = percent, fill = term)) +
  geom_col() +
  facet_wrap(~ term, ncol = 1) +
  scale_fill_manual(values = c("SBS6" = "#1F4E79", 
                               "SBS19" = "#6C757D",
                               "SBS8" = "#9E2A2B",
                               "TMB>=20" = "#000",
                               "SBS26" = "#4E6E58")) +
  ylab("% of tumour samples with exposure variable") +
  labs(title = 'B') +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = 'bold'),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22, face = 'bold'),
    axis.title.y = element_text(size = 22, face = 'bold'),
    strip.text = element_text(size = 22, face = 'bold'),
    axis.title.x = element_blank(),
    legend.position = "none")

png('summary_TMB_sig.png', width = 800, height = 1000)
p2
dev.off()


# rank in order of median TMB
tmb_median_ranks <- fread('TMB_summarystats.tsv') %>%
  arrange(desc(q3)) %>% filter(tumour != 'PANCANCER' & tumour != 'OTHER') %>%
  pull(tumour) %>% c('PANCANCER',.) %>% factor(., levels = unique(.))

# cancer_factor 
# or cancer sample size
cancer_factor <- fread("cohortsummary_bystatus.tsv") %>%
  # pull just the rows relating to cancer
  dplyr::slice(17:32) %>% 
  # remove the category 'OTHER' since this wasn't included in the tumour-stratified analysis
  filter(characteristic != 'OTHER') %>%
  # arrange in order of column 2 (frequencies)
  arrange(desc(`overall, n = N`)) %>% 
  # convert to a factor keeping the order they appear in this table as
  pull(characteristic) %>% c('PANCANCER',.) %>%
  factor(., levels = unique(.))

# Combine high TMB with signatures
# fully adjusted sigs
sigs_stratified <-fread('cox_signatures_tumour_stratified.tsv') %>%
  mutate(sensitivity_analysis = 'tumour-stratified') %>%
  mutate(stratum = cancer) %>%
  # remove results for any cancers where n_sig < 5 as regression model likely too unstable
  left_join(., freqs2, by = c('cancer', 'term')) %>% 
  mutate(across(c(HR, u95ci, l95ci, pval), ~case_when(n_sig < 5 ~ NA_real_,
                                                      .default = .))) %>%
  # now remove unnecessary columms
  select(-c(n_sig, percent_with_sig))

tmb_stratified <- fread('cox_TMB_tumour_stratified.tsv') %>%
  # just select tmb_binary_thresh20
  filter(term == 'tmb_binary_thresh20tmb_above_20') %>%
  mutate(term = 'TMB>=20') %>%
  # exclude any tumours where n TMB>20 is less than 5. I don't have these numbers outside the RE yet
  mutate(across(c(HR, l95ci, u95ci, pval), 
                ~case_when(cancer == 'HEPATOPANCREATOBILIARY'|
                             cancer == 'RENAL'|
                             cancer == 'SARCOMA'|
                             cancer == 'PROSTATE'|
                             cancer == 'HEAD_NECK'|
                             cancer == 'HAEM_ONC' ~ NA_real_, 
                           .default = .)))

#pancancer resutls
tmb_pancan <- fread('cox_TMB_fulladj.tsv') %>% 
  mutate(cancer = 'PANCANCER') %>%
  filter(term == 'tmb_binary_thresh20tmb_above_20') %>%
  mutate(term = 'TMB>=20')
sigs_pancan <- fread('cox_signatures_fulladj.tsv')  %>% 
  mutate(cancer = 'PANCANCER') %>% 
  mutate(term = gsub('signature_', 'SBS', term)) %>%
  filter(term %in% c('SBS6', 'SBS19', 'SBS26', 'SBS8'))

terms_to_plot <- c('TMB>=20', 'SBS6', 'SBS19', 'SBS26', 'SBS8')
df_plot1 <- rbind(tmb_pancan, sigs_pancan, 
                  sigs_stratified, 
                  tmb_stratified,
                  fill =T) %>%
  mutate(term = gsub('signature_', 'SBS', term)) %>%
  mutate(cancer = factor(cancer, levels = cancer_factor),
         # drop unused factor levels otherwise they still appear on graph
         term = factor(term, levels = terms_to_plot)) %>%
  # Arrange by cancer and gene for consistent ordering
  arrange(cancer, term) #%>%
  # remove cancers which have no results at all for any of the TMB/sigs
 # group_by(cancer)  %>% filter(!all(is.na(HR))) %>% ungroup()

p1 <- ggplot(df_plot1, aes(x = cancer, y = HR, color = term)) +
  geom_point(size = 3) +
  scale_y_continuous(trans = "log2") +
  geom_errorbar(aes(ymin = l95ci, ymax = u95ci), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ term, ncol = 1) +
  scale_color_manual(values = c("SBS6" = "#1F4E79", 
                                "SBS19" = "#6C757D",
                                "SBS8" = "#9E2A2B",
                                "TMB>=20" = "#000",
                                "SBS26" = "#4E6E58")) +
  ylab("Hazard Ratio for VTE") +
  labs(title = 'A') +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = 'bold'),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22, face = 'bold'),
    axis.title.y = element_text(size = 22, face = 'bold'),
    strip.text = element_text(size = 22, face = 'bold'),
    axis.title.x = element_blank(),
    legend.position = "none")
  

ggsave("./TMB_sigs_tumour_strat_forest.jpg", plot = p, width = 10, height = 12, dpi = 300)

# combined plot side by side
ggsave("./figure4.jpg", plot = p1 + p2, width = 20, height = 16, dpi = 300)

# Supp Table 10 germline PRS and thrombophilia SNPs #######
# Germline ####
prs_minadj <- fread('cox_PRS_minadj.tsv') %>% mutate(model = 'min_adj')
prs_fulladj <- fread('cox_PRS_fulladj.tsv') %>% mutate(model = 'full_adj') 
df1 <- fread('Cox_PRS_gene_interactions.tsv') %>% 
  mutate(model = 'gene:PRS interaction') %>% filter(str_detect(term, ':'))
df2 <- fread('cox_PRS_TMB_interactions.tsv') %>% 
  mutate(model = 'tmb:PRS interaction') %>% filter(str_detect(term, ':'))
df3 <- fread('cox_PRS_signature_interactions.tsv') %>% 
  mutate(model = 'SBS:PRS interaction') %>% filter(str_detect(term, ':'))
df4 <- fread('cox_germlineSNP_somatic_genes.tsv') %>%
  # filter just for the core results for rs6025 and rs1799963 and the interaction terms
  filter(term %in% c('rs6025_Het', 'rs6025_Hom', 'rs1799963_Het')|
           str_detect(term, ':')) %>%
  # change term labelling so it is clearer
  mutate(term = gsub('_Het', ' heterozygote', term),
         term = gsub('_Hom', ' homozygote', term),
         term = gsub('_genotype1', ' heterozygote', term))

df <-  rbind(prs_minadj, prs_fulladj, df1, df2, df3, df4, fill =T) %>%
  mutate(term = gsub('klarin_PRS_scaled', 'germline_PRS', term),
         term = gsub('klarin_PRS_scaled', 'germline_PRS', term),
         model_covars = gsub('klarin_PRS_scaled|prs', 'germline_PRS', model_covars),
         term = gsub('signature_', 'SBS', term),
         model= gsub('interaction', 'interaction_term', model),
         model_covars = gsub('signature_', 'SBS', model_covars)) %>%
  relocate(model, .before = HR) %>%
  mutate(across(c(HR, u95ci, l95ci), ~format(round(., digits = 2), nsmall =2))) %>%
  mutate(pval = formatC(pval, format = 'e', digits =2)) %>%
  mutate('HR [95% CI]' = paste0(HR, ' [', l95ci, '-', u95ci, ']')) %>%
  dplyr::select(term, model, 'HR [95% CI]', pval,  model_covars)

write.csv(df, './germline_results_formatted.csv', row.names= F)  

# Figure 5 PRS/PCDH15 and FVL PCDH15 interaction ######
## FVL freqs
freqs <- fread('SNP_freqs_FVL_PTG20210A.tsv')

data1 <- fread('cuminc_plot_PRS_PCDH15_cmprsk.tsv')

data2 <- fread('cuminc_plot_FVL_PCDH15_cmprsk.tsv')

prs_pcdh15 <- data1 %>%
      mutate(
        # already inverted surv to create cuminc inside RE but need to do same with CIs
        cum_inc_upper = 1 - lower,
        cum_inc_lower = 1 - upper,
        risk_level = case_when(
          grepl("high PRS", strata) ~ "high germline polygenic risk score",
          grepl("low PRS", strata) ~ "low germline polygenic risk score"),
        has_mutation = ifelse(grepl("PCDH15 mutated", strata), "PCDH15 mutated", 
                              "PCDH15 wildtype"))

fvl_pcdh15 <- data2 %>% mutate(
  cum_inc_upper = 1 - lower,
  cum_inc_lower = 1 - upper,
        risk_level = ifelse(grepl("FVL carrier", strata), "factor V Leiden heterozygote", 
                            "factor V Leiden wildtype"),
        has_mutation = ifelse(grepl("PCDH15 mutated", strata), "PCDH15 mutated", 
                              "PCDH15 wildtype")
      )
  
  
  # Plot
create_plot <- function(df,  title) {
  p <- ggplot(df %>% 
                # convert surv to a %
                mutate(surv = surv*100), aes(x = time, y = surv, color = risk_level, linetype = has_mutation,
                      fill = risk_level)) +
    geom_line(linewidth = 1) +
    scale_x_continuous(breaks = seq(0, max(df$time), 12)) +
    geom_vline(xintercept = 6, color = "gray60", linetype = "dotted", linewidth = 0.8) +
    scale_color_manual(values = c("low germline polygenic risk score" = "#1F4E79", "high germline polygenic risk score" = "#9E2A2B",
                                  "factor V Leiden heterozygote" = "#9E2A2B", "factor V Leiden wildtype" = "#1F4E79")) +
    scale_fill_manual(values = c("low germline polygenic risk score" = "#1F4E79", "high germline polygenic risk score" = "#9E2A2B",
                                 "factor V Leiden heterozygote" = "#9E2A2B", "factor V Leiden wildtype" = "#1F4E79")) +
    scale_linetype_manual(values = c("PCDH15 mutated" = "solid", 
                                     "PCDH15 wildtype" = "dashed"),
                          labels = c("solid: PCDH15 mutated",
                                     "dashed: PCDH15 wildtype")) +
    labs(title = title, x = "Time (months)", y = "VTE cumulative incidence (%)") +
    scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, 5),
                       minor_breaks = seq(0, 35, 1)) +
    labs(linetype = NULL, fill = NULL, colour = NULL) +
    guides(
      linetype = guide_legend(nrow = 2, byrow = TRUE, direction = "vertical"),
      fill = guide_legend(nrow = 2, byrow = TRUE, direction = "vertical"),
      color = guide_legend(nrow = 2, byrow = TRUE, direction = "vertical"),
    ) +
    theme_minimal() +
    theme(text = element_text(size = 20),
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
      legend.position = "top",
      legend.box = "vertical"
    )
  
  return(p)
}

# Usage
p1 <- create_plot(prs_pcdh15, 'A')
p2 <- create_plot(fvl_pcdh15, 'B')

ggsave("germline_somatic_interactions.jpeg", p1 + p2, 
       width = 210, height = 180, dpi = 300, units = "mm")

# cuminc at 6months
cuminc <- fvl_pcdh15 %>% filter(time %in% c(6, 12)) %>% 
  select(strata, time, surv, cum_inc_upper, cum_inc_lower) %>%
  mutate(across(c(surv, cum_inc_upper, cum_inc_lower), ~format(round(., digits=2), nsmall=2))) %>%
  pivot_wider(names_from = 'time', values_from = c(surv, cum_inc_upper, cum_inc_lower)) %>%
  mutate(cum.inc_6m = paste(surv_6, '[', cum_inc_lower_6, '-', cum_inc_upper_6, ']'),
         cum.inc_12m = paste(surv_12, '[', cum_inc_lower_12, '-', cum_inc_upper_12, ']'),
         stratum = strata) %>%
  select(stratum, starts_with('cum.inc')) %>%
  mutate(stratum = gsub('carrier', 'heterozygote', stratum),
         stratum = gsub('FVL', 'rs6025', stratum))

fwrite(cuminc, 'cuminc_6_12M_pcdh15_fvl.csv', row.names = F)  

# new figures added after peer review
# Figure 6 cuminc cancertype + mutation ######
# cuminc plots
filepaths <- list.files(".", 
                         pattern = "^cuminc_plot_Li.*\\.tsv$", full.names = TRUE)

# Read and combine all files
combined_data <- map_dfr(filepaths, function(file) {
  gene_name <- str_extract(file, "(?<=cuminc_plot_Li_cancer_category_).*(?=_cmprsk\\.tsv)")
  
  fread(file) %>%
    mutate(
      cuminc = surv, cuminc_ucl = upper, cuminc_lcl = lower,
      gene = gene_name,
           strata = str_remove_all(strata, "_") %>%
             str_replace("cancerplussomaticmutation", " +sm") %>%
             str_replace("cancernomutation", ""))
}) %>% select(gene, strata, cuminc, cuminc_lcl, cuminc_ucl, time) 
rm(filepaths)
  # remove the rows which just relate to cancer stage
ca_only <- combined_data %>% filter(is.na(gene)) %>%
  mutate(strata = factor(strata, 
      levels = c("Low risk",
                 "Intermediate risk", "High risk",
                 "Very high risk")))
genes_and_ca <- combined_data %>% filter(!is.na(gene)) %>%
  mutate(strata = factor(strata, 
                         levels = c("veryhighrisk +sm", "veryhighrisk",
                                    "highrisk +sm", "highrisk", "intermediaterisk +sm",
                                    "intermediaterisk", "lowrisk +sm", "lowrisk")))

cindex <- fread('cindex_li_cancer_category_plus_somatic_mutation.txt') %>%
  dplyr::rename(gene = V1) %>% 
  mutate(label = paste0(gene, " + cancer-risk group \n (C-index: ", round(cindex, 2), ")"))

risk_colors <- c(
  "veryhighrisk +sm" = "#8B0000",      # dark red
  "veryhighrisk" = "#CD5C5C",   # light red
  "highrisk +sm" = "#FF8C00",       # dark orange
  "highrisk" = "#FFB84D",    # light orange
  "intermediaterisk +sm" = "#1E90FF",       # dark blue
  "intermediaterisk" = "#87CEEB",    # light blue
  "lowrisk +sm" = "#228B22",       # dark green
  "lowrisk" = "#90EE90"     # light green
)
p <- ggplot(genes_and_ca %>% left_join(., cindex, by = 'gene') %>%
              # get 6 month cuminc
              filter(time == 6), 
       aes(x = strata, y = cuminc*100, ymin = cuminc_lcl*100,
           ymax = cuminc_ucl*100,
           colour = strata)) +
  geom_pointrange() +
  scale_colour_manual(values = risk_colors) +
  scale_y_continuous(breaks = seq(0, 12.5, 2.5)) +
  facet_wrap(~label) +
  labs(x = "", y = "6-month VTE cumulative incidence (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 14, face = 'bold'),
        axis.title.y = element_text(size = 18, face = 'bold'),
        strip.text = element_text(size = 16, face = 'bold'),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA))


ggsave("cuminc_6M_cancer_category_gene.jpeg", p, 
       width = 190, height = 210, dpi = 300, units = "mm")

# Supp Table 11 cuminc stratified by tumour type and stage ##########
# number of participants in each category
filepaths <- list.files(".", 
                        pattern = "^freqs_li.*\\.tsv$", full.names = TRUE)

# Read and combine all files
combined_freqs <- map_dfr(filepaths, function(file) {
  gene_name <- str_extract(file, "(?<=freqs_li_cancer_risk_grp_plus).*(?=\\.tsv)")
  fread(file) %>% mutate(gene = gene_name)
  }) %>%
  mutate(model = case_when(is.na(gene) ~ 'cancer risk category',
                           !is.na(gene) ~ paste('cancer risk category +', gene))) %>%
  dplyr::rename(stratum = Var1) %>%
  mutate(stratum = str_remove_all(stratum, "_") %>%
  str_replace("cancerplussomaticmutation", " +sm") %>%
  str_replace("cancernomutation", ""))

cuminc_6M <- rbind(ca_only, genes_and_ca) %>%
  filter(time %in% c(6,12)) %>%
  mutate(model = case_when(is.na(gene) ~ 'cancer risk category',
                           !is.na(gene) ~ paste('cancer risk category +', gene))) %>%
pivot_wider(names_from = 'time', values_from = c(cuminc, cuminc_lcl, cuminc_ucl)) %>%
  mutate(across(starts_with('cuminc'), ~format(round(., digits =3), nsmall =2)),
        cum.inc_6m = paste(cuminc_6, '[', cuminc_lcl_6, '-', cuminc_ucl_6, ']'),
         cum.inc_12m = paste(cuminc_12, '[', cuminc_lcl_12, '-', cuminc_ucl_12, ']'),
         stratum = strata) %>%
  # join with the frequencies table
  left_join(., combined_freqs, by = c('model', 'stratum')) %>%
  # N = total cohort size
  mutate(`freq (%)` = paste0(Freq, ' (', round((Freq/N)*100, digits =1), '%)'),
         stratum = gsub('Low risk', 'Low risk (all other tumour types)', stratum),
         stratum = gsub('Intermediate risk', 'Intermediate risk (colorectal/intestinal tumours)', stratum),
         stratum = gsub('High risk', 'High risk (brain, lung, gynaecological, genitourinary, sarcoma and high-grade haematological tumours)', stratum),
         stratum = gsub('Very high risk', 'Very high risk pancreatobiliary and gastro-oesophageal tumours', stratum),
         stratum = gsub('veryhighrisk', 'Very high risk', stratum),
         stratum = gsub('highrisk', 'High risk', stratum),
         stratum = gsub('intermediaterisk', 'Intermediate risk', stratum),
         stratum = gsub('lowrisk', 'Low risk', stratum),
         ) %>%
  select(model, stratum, `freq (%)`, starts_with('cum.inc'))
  
# do a wide format version for the gene + cancer types
cuminc_6M_wide <- cuminc_6M %>% 
  # filter(str_detect(model, '+')) %>%
  pivot_wider(names_from = 'model', values_from = c(`freq (%)`, starts_with('cum.inc'))) %>%
  select(stratum, `freq (%)_cancer risk category`, 
         `cum.inc_6m_cancer risk category`, `cum.inc_12m_cancer risk category`,
         contains('CDKN2A'), contains('KRAS'), contains('PCDH15'), contains('TP53'))

fwrite(cuminc_6M_wide, 'cuminc_6M_12M_cancer_type_mutation.csv', row.names =F)  


## cuminc cancerstage + mutation ######
# cuminc plots
filepaths <- list.files(".", 
                        pattern = "^cuminc_plot_stage.*\\.tsv$", full.names = TRUE)

# Read and combine all files
combined_data <- map_dfr(filepaths, function(file) {
  gene_name <- str_extract(file, "(?<=cuminc_plot_stage_).*(?=_cmprsk\\.tsv)")
  
  fread(file) %>%
    mutate(
      cuminc = surv, cuminc_ucl = upper, cuminc_lcl = lower,gene = gene_name,
           strata = str_remove_all(strata, "_") %>%
             str_replace("cancerplussomaticmutation", " +sm") %>%
             str_replace("cancernomutation", ""))
}) %>% select(gene, strata, cuminc, cuminc_lcl, cuminc_ucl, time) %>%
  mutate(strata = gsub('earlystagecancer|earlystage|early', 'early stage', strata),
         strata = gsub('latestagecancer|latestage|late', 'late stage', strata))
rm(filepaths)
# remove the rows which just relate to cancer stage
genes_and_stage <- combined_data %>% 
  mutate(strata = factor(strata, 
                         levels = c("late stage +sm",
                                    "late stage",
                                    "early stage +sm",
                                    "early stage"
                                   )))

cindex <- fread('cindex_stage_plus_somatic_mutation.tsv') %>%
  dplyr::rename(gene = V1) %>% 
  mutate(label = paste0(gene, " + stage \n (C-index: ", round(cindex, 2), ")"))

risk_colors <- c(
  "late stage +sm" = "#8B0000",      # dark red
  "late stage" = "#CD5C5C",   # light red
 "early stage +sm" = "#1E90FF",       # dark blue
  "early stage" = "#87CEEB"    # light blue
)

# Supp Figure 7 cuminc 6months by stage and gene ########
p <- ggplot(genes_and_stage %>% left_join(., cindex, by = 'gene') %>%
              # remove stage only
              filter(!is.na(gene)|gene != 'only') %>%
              # get 6 month cuminc
              filter(time == 6), 
            aes(x = strata, y = cuminc*100, ymin = cuminc_lcl*100,
                ymax = cuminc_ucl*100,
                colour = strata)) +
  geom_pointrange() +
  scale_colour_manual(values = risk_colors) +
  scale_y_continuous(breaks = seq(0, 12.5, 2.5)) +
  facet_wrap(~label) +
  labs(x = "", y = "6-month VTE cumulative incidence (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 14, face = 'bold'),
        axis.title.y = element_text(size = 18, face = 'bold'),
        strip.text = element_text(size = 16, face = 'bold'),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA))

ggsave("cuminc_6M_cancer_stage_gene.jpeg", p, 
       width = 190, height = 210, dpi = 300, units = "mm")

# number of participants in each category
filepaths <- list.files(".", 
                        pattern = "^freqs_stage.*\\.tsv$", full.names = TRUE)

# Read and combine all files
combined_freqs <- map_dfr(filepaths, function(file) {
  gene_name <- str_extract(file, "(?<=freqs_stage_plus).*(?=\\.tsv)")
  fread(file) %>% mutate(gene = gene_name,
                         strata = Var1,
                         strata = str_remove_all(strata, "_") %>%
                         str_replace("cancerplussomaticmutation", " +sm") %>%
                         str_replace("cancernomutation", ""))
}) %>% select(gene, strata, Freq) %>%
  mutate(strata = gsub('earlystagecancer|earlystage|early', 'early stage', strata),
         strata = gsub('latestagecancer|latestage|late', 'late stage', strata),
         gene = case_when(is.na(gene) ~ 'only', .default =gene))

cuminc_6M <- genes_and_stage %>%
  filter(time %in% c(6,12)) %>%
  mutate(model = case_when(is.na(gene) ~ 'cancer stage',
                           !is.na(gene) ~ paste('cancer stage +', gene))) %>%
  pivot_wider(names_from = 'time', values_from = c(cuminc, cuminc_lcl, cuminc_ucl)) %>%
  mutate(across(starts_with('cuminc'), ~format(round(., digits =3), nsmall =2)),
         cum.inc_6m = paste(cuminc_6, '[', cuminc_lcl_6, '-', cuminc_ucl_6, ']'),
         cum.inc_12m = paste(cuminc_12, '[', cuminc_lcl_12, '-', cuminc_ucl_12, ']')
         ) %>%
  # join with the frequencies table
  left_join(., combined_freqs,
            by = c('gene', 'strata')) %>%
  # N = total cohort size
  mutate(`freq (%)` = paste0(Freq, ' (', round((Freq/N)*100, digits =1), '%)'),
         strata = gsub('earlystagecancer|earlystage|early', 'early', strata),
         strata = gsub('latestagecancer|latestage|late', 'late', strata),
         stratum = strata
  ) %>%
  select(model, stratum, strata, `freq (%)`, starts_with('cum.inc'))

# do a wide format version for the gene + cancer types
cuminc_6M_wide <- cuminc_6M %>% 
  # filter(str_detect(model, '+')) %>%
  mutate(model = gsub(' \\+ only', '', model)) %>%
  pivot_wider(names_from = 'model', 
              values_from = c(`freq (%)`, starts_with('cum.inc'))) %>%
  select(stratum, `freq (%)_cancer stage`, 
         `cum.inc_6m_cancer stage`, `cum.inc_12m_cancer stage`,
         contains('CDKN2A'), contains('KRAS'), contains('PCDH15'), contains('TP53'))

fwrite(cuminc_6M_wide, 'cuminc_6M_12M_cancer_stage_mutation.csv', row.names =F)  


# Tumour interaction models (main text) #####
anova_genes_tumour <- fread('anova_tumour_gene_interactions.tsv')
anova_sig_tumour <- fread('anova_tumour_signature_interactions.tsv')
anova_tmb_tumour <- fread('anova_tumour_TMB_interactions.tsv')
int_genes_tumour <- fread('interaction_terms_tumour_gene.tsv')

# Supp Figure 8 replication of previously published associations #######
genes_all <- fread('genes_all.csv')
genes_prev_reported_higher_VTE <- c('CDKN2B','CTNNB1','KEAP1',
  'MET','STK11','TP53','KRAS',
  'EGFR','PIK2CA','CDKN2A','PTEN','NF1',
  'ALK', 'ROS1')
genes_prev_reported_lower_VTE <- c('IDH1', 'SETD2')
# Filter gene list
plot_data <- genes_all %>%
  arrange(gene)%>%
  pivot_longer(cols = c(HR_minadj, HR_fulladj),
               names_to = "model",
               values_to = "HR",
               names_prefix = "HR_") %>%
  mutate(l95ci = if_else(model == "minadj", l95ci_minadj, l95ci_fulladj),
         u95ci = if_else(model == "minadj", u95ci_minadj, u95ci_fulladj),
         model = factor(model, levels = c("minadj", "fulladj"),
                        labels = c("min.adj", "full.adj"))) %>%
  mutate(model = factor(model, levels = c('full.adj', 'min.adj'))) %>%
  arrange(gene, model) 

# genes prev reported to increase VTE
p1 <- ggplot(plot_data %>%
              filter(gene %in% genes_prev_reported_higher_VTE), aes(y = model, x = HR, color = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = l95ci, xmax = u95ci), height = 0.2) +
  coord_trans(x = "log2") +
  geom_text(aes(label = sprintf("HR: %.2f (%.2f-%.2f)", 
                                HR, l95ci, u95ci)),
            x = Inf, hjust = 1, size = 3, color = "black") +
  facet_wrap(~gene, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_x_continuous(
    breaks = c( 0.25, 0.5, 1, 2, 4),
    labels = c("0.25", "0.5", "1", "2", "4"),
    limits = c(0.25,4)) +
  scale_color_manual(values = c("#1F4E79", "black")) +
  labs(title = 'A',
       x = "Hazard Ratio", y = NULL) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.text.y.left = element_text(size = 16, face = 'bold', angle = 0, hjust = 1),
        strip.background = element_rect(fill = "grey90"),
        panel.grid = element_blank())

ggsave("prev_reports_higher.jpeg", p1, 
       width = 190, height = 210, dpi = 300, units = "mm")

# genes prev reported to decrease VTE
p2 <- ggplot(plot_data %>% 
              filter(gene %in% genes_prev_reported_lower_VTE), aes(y = model, x = HR, color = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = l95ci, xmax = u95ci), height = 0.2) +
  coord_trans(x = "log2") +
  geom_text(aes(label = sprintf("HR: %.2f (%.2f-%.2f)", 
                                HR, l95ci, u95ci)),
            x = Inf, hjust = 1, size = 3, color = "black") +
  facet_wrap(~gene, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_x_continuous(
    breaks = c( 0.25, 0.5, 1, 2, 4),
    labels = c("0.25", "0.5", "1", "2", "4"),
    limits = c(0.25,4)) +
  scale_color_manual(values = c("#1F4E79", "black")) +
  labs(title = 'B',
       x = "Hazard Ratio", y = NULL) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.y.left = element_text(size = 16, face = 'bold', angle = 0, hjust = 1),
        strip.background = element_rect(fill = "grey90"),
        panel.grid = element_blank())

p3 <- (p1 / p2) + plot_layout(widths = 1, heights = c(12, 2))

ggsave("prev_reported_genes.jpeg", p3, 
       width = 190, height = 210, dpi = 300, units = "mm")
