# TMT-SILAC analysis script. Call and run this at the start of the dashboard.
# ==============================================================================
# 
# Libraries
library(tidyverse)
library(readr)
library(robustbase)
# Also requires plyr, DescTools, matrixStats, ggcorrplot

########## IMPORTANT PARAMETERS - decide were to call them ########################
# MOST USUAL: If running this script from the dashboard these are defined in that markdown file

# If sourcing as a background job directly in RStudio you can just set them here
# working_dir <- "toy_data"
# Times <- c(21,16,12,4,2,1,0.25,0) * 24

# If just calling from the command line, you can pass them as arguments using this
# args = commandArgs(trailingOnly=TRUE)
# 
# if (length(args) < 2) {
#   stop("Please supply 2 arguments: 1. data directory, 2. timepoints (e.g. '21,16,12,4,2,1,0.25,0')", call.=FALSE)
# } else if (length(args) == 2) {
#   # default output file
#   working_dir <- args[1]
#   Times <- as.numeric(unlist(strsplit(args[2],",")))
# }
##################################################################################

# Import TMT and SILAC tables ----

source("R_functions_revised/importIP2.R")
source("R_functions_revised/helper_functions/readTMT.R")
source("R_functions_revised/helper_functions/readSILAC.R")

# Create design table, alternatively import it.
files <- list.files(path = paste0(working_dir,"/IP2_data"), full.names = T, recursive = T)

samples <- t(as.data.frame(str_split(str_match(gsub(pattern = "^.*?/IP2_data/",replacement = "",x = files),pattern = "(.*?).(csv|xls)")[,2],pattern = "_")))

design <- as.data.frame(cbind(files = files, localization = samples[,1], datatype = tolower(samples[,2])),stringsAsFactors = F, row.names = seq.int(1:length(files)))

outputdir = paste0(working_dir, "/output_dir")
dir.create(outputdir)

raw_data <- importIP2(design$files, design$datatype, design$localization, outputdir)

rm(files, samples, design)

# Remove empty TMT channels and/or rename columns according to timecourse
source("R_functions_revised/rmTMTtagIP2.R")
raw_data <- rmTMTtagIP2(raw_data, Times)


# Filter TMT and SILAC tables ----

# Filter out common contaminants and poor data (blank timepoints, etc)
source("R_functions_revised/filterIP2.R")
source("R_functions_revised/helper_functions/filterTMT.R")
source("R_functions_revised/helper_functions/filterSILACv2.R")

filtered_data <- filterIP2(data = raw_data, filter_description = c("Sus scrofa", "Staphylococcus aureus"), outputdir = outputdir)

# Keep track of proteins present:
# 
# Number of proteins after 1st filter (not counted before this to avoid contaminants and reverse peps):

# Total
filtered_prot_num <- filtered_data[["SILAC"]] %>% ungroup() %>% select(uniprot) %>%
    unique() %>% count() %>% as.numeric(.)

# Per cell compartment
Prot_CC <- filtered_data[["SILAC"]] %>% ungroup() %>% 
    select(uniprot, Localization) %>% 
    unique() %>% group_by(Localization) %>% count()

# Create table to track info per pipeline step:
pip_step_prot <- setNames(data.frame(matrix(ncol = length(Prot_CC$Localization)+2, nrow = 0)), c("pipeline_step",Prot_CC$Localization, "total"))

# Add data:
pip_step_prot[1,] <- c("starting filtered proteins", Prot_CC$n, filtered_prot_num)
rm(Prot_CC, filtered_prot_num)

# Perform row and sum normalization of the data
source("R_functions_revised/normIP2_v2.R")

not_normalized_tmt <- normIP2(filtered_data, row_norm = FALSE, sum_norm = FALSE)
all_normalized_tmt <- normIP2(data = filtered_data, row_norm = TRUE, sum_norm = TRUE, outputdir = outputdir)

# Save barplot before and plot after normalization

source("R_functions_revised/plotting_functions.R")
p_no_norm <- pcnt_label_plot(normalized_tmt = not_normalized_tmt, F)
p_all_norm <- pcnt_label_plot(normalized_tmt = all_normalized_tmt, F)
rm(not_normalized_tmt)
# CC barplots only after normalization:
p_all_norm_CC <- pcnt_label_plot(normalized_tmt = all_normalized_tmt, T)

# Number of proteins after normalization:
# Total
norm_prot_num <- as.numeric(all_normalized_tmt %>% ungroup() %>% select(uniprot) %>% unique() %>% count())
# Per cell compartment:
Prot_CC <- all_normalized_tmt %>% ungroup() %>% select(uniprot, Localization) %>% unique() %>% group_by(Localization) %>% count()
# Add to table:
pip_step_prot <- rbind(pip_step_prot, c("normalized proteins", Prot_CC$n, norm_prot_num))
rm(Prot_CC,norm_prot_num)

# Calculate fold-change for each petide across timepoints:
source("R_functions_revised/FC_heavy_light.R")
all_norm_fc <- FC_calc(all_normalized_tmt,Times)

# Number of light proteins to be used for decay rate calculation:
lightfc_prot_num <- as.numeric(all_norm_fc %>% ungroup() %>% filter(Type == "TMT_light") %>% select(uniprot) %>% unique() %>% count())
Prot_CC <- all_norm_fc %>% ungroup() %>% filter(Type == "TMT_light") %>% select(uniprot,Localization) %>% unique() %>% group_by(Localization) %>% count()
pip_step_prot <- rbind(pip_step_prot, c("light proteins (after FC)", Prot_CC$n, lightfc_prot_num))
rm(Prot_CC, lightfc_prot_num)
# Number of heavy proteins to be used for synthesis rate calculation
heavyfc_prot_num <- as.numeric(all_norm_fc %>% ungroup() %>% filter(Type == "TMT_heavy") %>% select(uniprot) %>% unique() %>% count())
Prot_CC <- all_norm_fc %>% ungroup() %>% filter(Type == "TMT_heavy") %>% select(uniprot, Localization) %>% unique() %>% group_by(Localization) %>% count()
pip_step_prot <- rbind(pip_step_prot, c("heavy proteins (after FC)", Prot_CC$n, heavyfc_prot_num))
rm(Prot_CC, heavyfc_prot_num)

# Save lineplot of synthesis and decay overall
############################ ERROR HERE! ###########################
source("R_functions_revised/plotting_functions.R")
p_syndec_all <- proteome_mean_plot(all_norm_fc,F)
# Save line plot of synthesis and decay for each cell compartment
p_syndec_CC <- proteome_mean_plot(all_norm_fc,T)

# You can also have plots with the NLS fit instead of lines:
p_syndec_all_nls <- proteome_mean_plot2(all_norm_fc,F)
p_syndec_CC_nls <- proteome_mean_plot2(all_norm_fc,T)
#####################################################################
# Guess Thalf ----

# Guestimate half-life based on overall percent heavy:
# Using the regular exponential decay curve. I gave up on the monte carlo simulation here. I'n not too sure how that can be useful since averaging those values will simple tend to the actual direct estimate.

source("R_functions_revised/guess_Thalf_pctHeavy_fixed.R")
all_norm_fc <- guess_Thalf_pctHeavy(Times, all_norm_fc)

# Not sure it makes sense to also count the number of proteins here... Probably not since they will all get an estimate. Maybe I could count what ended up outside the min and max values
# IMPORTANT: Given the assumption of a time point being 100% light always, the max Pct Heavy will determine the minimal Thalf that could be estimated and this will mainly depend on the number of timepoints present. This could be why we can't see the short lived protein bump we see in the nls fits.
# This also does not match what we see experimentally because some of the peptides have higher pect heavy than the theoretically allowed amount. We also see a bit of heavy in the time zero. Not sure why.

rm(raw_data, filtered_data,all_normalized_tmt)

# Fit decay and synthesis models, peptide level:
# Fit port model for synthesis and decay. Use the regular fit and the robust version.
# No outliers are being aliminated right now.

# Fit NLS models ----

source("R_functions_revised/nls_mods.R")
groupings <- c("Unique", "Sequence", "uniprot", "description", "prot_name", "gene_sym", "Localization", "Type")
# Careful with the savemod option, this takes a long time (over 1h) and models are huge, over 5GB at least. I have a feeling saving might actually take longer than re-running it.
All_mods <- nls_mods(all_norm_fc = all_norm_fc, groupings = groupings, savemod = F, outputdir = outputdir, Times = Times)

# Separate data for plotting from model results:
pred_data <- All_mods$NLS_preds
All_mods <- All_mods$NLS_results
rm(all_norm_fc)

# Replace poor fits with NA, then filter data, I think I'll use 0.8 correlation and 0.25 rmse for now
All_mods <- All_mods %>% 
    mutate(alpha_port = ifelse(cor_port < 0.8 | rmse_port > 0.25, NA, alpha_port),
           alpha_robMport = ifelse(cor_robMport < 0.8 | rmse_robMport > 0.25, NA, alpha_robMport)) %>%
    rename(alpha_guess = guess_alpha) %>% select(one_of(groupings), starts_with("alpha_")) %>% 
    pivot_longer(cols = starts_with("alpha"), names_to = "model",names_prefix = "alpha_",values_to = "alpha") %>% 
    filter(!is.na(alpha)) %>% mutate(Thalf = log(2)/alpha)

write_csv(All_mods,paste0(outputdir,"/NLS_models_results_filter.csv"))

# Number of proteins with fits per model:
filt_fit <- All_mods %>% select(uniprot, Type, model) %>% unique() %>% group_by(Type, model) %>% count()

# Number of proteins with fits per model per CC:
filt_fit_CC <- All_mods %>% select(uniprot, Type, Localization, model) %>% unique() %>% group_by(Localization,Type, model) %>% count()
filt_fit_CC <- filt_fit_CC %>% unite("pipeline_step",Type:model) %>% mutate(pipeline_step = paste0("filtered_fits_", gsub("TMT_","",pipeline_step))) %>% 
    pivot_wider(names_from = Localization, values_from = n)
filt_fit_CC$total <- filt_fit$n

pip_step_prot <- rbind(pip_step_prot, filt_fit_CC)
rm(filt_fit, filt_fit_CC)

# Aggregate per protein ----

### Calculate the median, mean, sd and IQR, for the Thalf. Add 
### peptide number count for each protein, for each model, per cell compartment. 
### Filter only unique peptides. Filter proteins with < 2 peptides.
### Also filter Inf Thalf peptides? Later...

All_mods_U <- All_mods %>% filter(Unique == T)

p_funs <- list(pep_num = ~ n(), alpha_med = ~ median(., na.rm = TRUE), alpha_iqr = ~ IQR(., na.rm = TRUE),alpha_avg = ~ mean(., na.rm = TRUE), alpha_sd = ~ sd(., na.rm = TRUE))

prot_alpha_CC_U <- All_mods_U %>% select(-Thalf) %>% 
    group_by_at(c("uniprot", "description", "prot_name", "gene_sym", "Localization","Type","model")) %>% 
    summarise_at(.vars = vars("alpha"), p_funs)

prot_alpha_CC <- All_mods %>% select(-Thalf) %>% 
    group_by_at(c("uniprot", "description", "prot_name", "gene_sym", "Localization","Type","model")) %>% 
    summarise_at(.vars = vars("alpha"), p_funs)

prot_alpha_CC_U_2 <- prot_alpha_CC_U %>% filter(pep_num >1) 
prot_alpha_CC_2 <- prot_alpha_CC %>% filter(pep_num >1)

# Count number of proteins left per model, per syn/dec per cell compartment only U or all.

filt_fit_U <- prot_alpha_CC_U_2 %>% ungroup() %>% select(uniprot, Type, model) %>% unique() %>% group_by(Type, model) %>% count()
filt_fit <- prot_alpha_CC_2 %>% ungroup() %>% select(uniprot, Type, model) %>% unique() %>% group_by(Type, model) %>% count()

# Number of proteins with fits per model per CC:

filt_fit_CC <- prot_alpha_CC_2 %>% ungroup() %>% select(uniprot, Type, Localization, model) %>% unique() %>% group_by(Localization,Type, model) %>% count()
filt_fit_CC <- filt_fit_CC %>% unite("pipeline_step",Type:model) %>% mutate(pipeline_step = paste0("prot_2peps_", gsub("TMT_","",pipeline_step))) %>% 
    pivot_wider(names_from = Localization, values_from = n)
filt_fit_CC$total <- filt_fit$n

filt_fit_CC_U <- prot_alpha_CC_U_2 %>% ungroup() %>% select(uniprot, Type, Localization, model) %>% unique() %>% group_by(Localization,Type, model) %>% count()
filt_fit_CC_U <- filt_fit_CC_U %>% unite("pipeline_step",Type:model) %>% mutate(pipeline_step = paste0("prot_U_2peps_", gsub("TMT_","",pipeline_step))) %>% 
    pivot_wider(names_from = Localization, values_from = n)
filt_fit_CC_U$total <- filt_fit_U$n

pip_step_prot <- rbind(pip_step_prot, filt_fit_CC, filt_fit_CC_U)
rm(filt_fit, filt_fit_CC, filt_fit_CC_U, filt_fit_U)

### Calculate the median. mean, sd, IQR and peptide number for each protein, for 
### each model, regardless of cell compartment. Filter unique peptides and 
### proteins with > 2 peptides.
prot_alpha_U <- All_mods_U %>% select(-Thalf) %>% 
    group_by_at(c("uniprot", "description", "prot_name", "gene_sym", "Type","model")) %>% 
    summarise_at(.vars = vars("alpha"), p_funs)

prot_alpha <- All_mods %>% select(-Thalf) %>% 
    group_by_at(c("uniprot", "description", "prot_name", "gene_sym", "Type","model")) %>% 
    summarise_at(.vars = vars("alpha"), p_funs)

prot_alpha_U_2 <- prot_alpha_U %>% filter(pep_num >1) 
prot_alpha_2 <- prot_alpha %>% filter(pep_num >1)

################# Using weighted median to aggregate data with bigger weight on unique peptides 

count_peps_CC <- All_mods %>% 
    group_by_at(c("uniprot", "description", "prot_name", "gene_sym", "Localization","Type","model")) %>%
    count(uniprot, Unique) %>% 
    mutate(unique_pep = if_else(Unique == TRUE,"U","not_U"))

count_peps_CC <- count_peps_CC %>% select(-Unique) %>% 
    pivot_wider(names_from = unique_pep, names_prefix = "pep_num_",values_from = n, values_fill = list(n = 0))

count_peps <- All_mods %>% 
    group_by_at(c("uniprot", "description", "prot_name", "gene_sym","Type","model")) %>%
    count(uniprot, Unique) %>% 
    mutate(unique_pep = if_else(Unique == TRUE,"U","not_U"))

count_peps <- count_peps %>% select(-Unique) %>% 
    pivot_wider(names_from = unique_pep, names_prefix = "pep_num_",values_from = n, values_fill = list(n = 0))

# Create a new column giving a weight of 10 to unique peptides and of 1 to non-uniques. Calculate the weighted median. Join this with the unique and non-unique peptide counts. Lastly join this with the other estimates, sd, iqr, etc.

All_mods_2 <- All_mods %>% 
    mutate(md_wts = if_else(Unique == T,10,1)) %>% 
    group_by_at(c("uniprot", "description", "prot_name", "gene_sym", "Type","model")) %>% 
    summarise(alpha_wt_med = matrixStats::weightedMedian(x = alpha, w = md_wts, na.rm = TRUE))
All_mods_2_CC <- All_mods %>% 
    mutate(md_wts = if_else(Unique == T,10,1)) %>% 
    group_by_at(c("uniprot", "description", "prot_name", "gene_sym", "Type","model","Localization")) %>% 
    summarise(alpha_wt_med = matrixStats::weightedMedian(x = alpha, w = md_wts, na.rm = TRUE))

All_mods_2 <- full_join(All_mods_2, count_peps, by = c("uniprot", "description", "prot_name", "gene_sym", "Type","model"))
All_mods_2_CC <- full_join(All_mods_2_CC, count_peps_CC, by = c("uniprot", "description", "prot_name", "gene_sym", "Type","model","Localization"))

All_mods_2 <- full_join(All_mods_2, prot_alpha, by = c("uniprot", "description", "prot_name", "gene_sym", "Type","model")) 
All_mods_2 <- All_mods_2 %>% mutate(Thalf_wt_med = log(2)/alpha_wt_med)
All_mods_2_CC <- full_join(All_mods_2_CC, prot_alpha_CC, by = c("uniprot", "description", "prot_name", "gene_sym", "Type","model","Localization"))
All_mods_2_CC <- All_mods_2_CC %>% mutate(Thalf_wt_med = log(2)/alpha_wt_med)

################################## IMPORTANT OUTPUT #########################################
# SAVE THE All_mods_2 AND All_mods_2_CC TO USE LATER IN THE COMBINED ANALYSIS BETWEEN SPECIES

save(All_mods_2_CC, file = paste0(outputdir,"/CC_protein_alpha_data.RData"))
save(All_mods_2, file = paste0(outputdir, "/protein_alpha_data.RData"))
#############################################################################################

# Save plots ----

# Half-life estimate summaries and plot:
### Save histograms of half-life distributions for synthesis and decay per model:

p_thalf_hist <- All_mods_2 %>% 
    ggplot(aes(x = Thalf_wt_med)) + geom_histogram() + 
    scale_x_continuous(trans = "log2") + 
    xlab("Median protein half-life (hours)") +
    facet_grid(Type ~ model) +
    theme_bw()

p_thalf_dens <- All_mods_2 %>% 
    ggplot(aes(x = Thalf_wt_med)) + geom_density(fill = "grey") + 
    scale_x_continuous(trans = "log2") +
    xlab("Median protein half-life (hours)") +
    facet_grid(Type ~ model) +
    theme_bw()

p_thalf_box <- All_mods_2 %>% 
    ggplot(aes(x = model, y = Thalf_wt_med, fill = Type)) + geom_boxplot() + 
    scale_y_continuous(trans = "log2") +
    ylab("Median protein half-life (hours)") +
    theme_bw()

# Same per cell compartment:
p_thalf_hist_CC <- All_mods_2_CC %>% 
    ggplot(aes(x = Thalf_wt_med, fill = Type)) + 
    geom_histogram(position = "dodge") + 
    scale_x_continuous(trans = "log2") + 
    xlab("Median protein half-life (hours)") +
    facet_grid(model ~ Localization) +
    theme_bw() + 
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, hjust = 1))

p_thalf_dens_CC <- All_mods_2_CC %>% 
    ggplot(aes(x = Thalf_wt_med, color = Type, fill = Type)) + 
    geom_density(alpha=.5) + 
    scale_x_continuous(trans = "log2") +
    xlab("Median protein half-life (hours)") +
    facet_grid(model ~ Localization) +
    theme_bw() + 
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, hjust = 1))

p_thalf_box_CC <- All_mods_2_CC %>% 
    ggplot(aes(x = Localization, y = Thalf_wt_med, fill = Type)) + 
    geom_boxplot() + 
    scale_y_continuous(trans = "log2") +
    ylab("Median protein half-life (hours)") +
    facet_grid(. ~ model) +
    theme_bw() + 
    theme(legend.position = "bottom")

### Calculate mean, median, sd and IQR for all proteins together:
# all proteins:
p_funs <- list(prot_num_tt = ~ n(), Thalf_wt_med_all = ~ median(., na.rm = TRUE), Thalf_iqr_all = ~ IQR(., na.rm = TRUE), Thalf_avg_all = ~ mean(., na.rm = TRUE), Thalf_sd_all = ~ sd(., na.rm = TRUE))

prot_median_all <- All_mods_2 %>% ungroup() %>% 
    select(Type, model, Thalf_wt_med) %>% 
    group_by(Type, model) %>%
    summarise_at(.vars = vars(starts_with("Thalf")), p_funs)

prot_median_all_noInf <- All_mods_2 %>% ungroup() %>% 
    select(Type, model, Thalf_wt_med) %>% 
    filter(Thalf_wt_med != Inf) %>% 
    group_by(Type, model) %>%
    summarise_at(.vars = vars(starts_with("Thalf")), p_funs)

# per cellular compartment:
prot_median_all_CC <- All_mods_2_CC %>% ungroup() %>% 
    select(Localization, Type, model, Thalf_wt_med) %>% 
    group_by(Localization, Type, model) %>%
    summarise_at(.vars = vars(starts_with("Thalf")), p_funs)

prot_median_all_noInf_CC <- All_mods_2_CC %>% ungroup() %>% 
    select(Localization, Type, model, Thalf_wt_med) %>% 
    filter(Thalf_wt_med != Inf) %>% 
    group_by(Localization, Type, model) %>%
    summarise_at(.vars = vars(starts_with("Thalf")), p_funs)

rm(p_funs,count_peps, count_peps_CC, prot_alpha, prot_alpha_CC)

# How do the half-lives correlate between the synthesis and decay for the 3 different models?
p_syn_deg_thalf <- All_mods_2 %>%
    filter(Thalf_wt_med != Inf) %>% 
    select(Type, uniprot, Thalf_wt_med, model) %>% 
    pivot_wider(names_from = Type, values_from = Thalf_wt_med) %>% 
    ggplot(aes(log2(TMT_heavy), log2(TMT_light))) + geom_point() +
    xlim(0,15) + ylim(0,15) +
    geom_smooth(method = "lm") +
    facet_wrap(.~ model)
# Weird issue with 2 lines containing the same Uniprot ID but from different sources. I might need to filter all the SWISS-PROT stuff. This means I loose the GFP, which should be fine.
p_cor <- All_mods_2 %>% ungroup() %>% 
    filter(Thalf_wt_med != Inf) %>%
    filter(!grepl("SWISS-PROT",prot_name)) %>% 
    select(Type, uniprot, Thalf_wt_med, model) %>% 
    mutate(Type = if_else(Type == "TMT_light","deg","syn")) %>% 
    pivot_wider(id_cols = uniprot, names_from = c(model, Type), values_from = Thalf_wt_med) %>% 
    select(-uniprot) %>% cor(use = "pairwise.complete.obs", method = "pearson")

s_cor <- All_mods_2 %>% ungroup() %>% 
    filter(Thalf_wt_med != Inf) %>% 
    filter(!grepl("SWISS-PROT",prot_name)) %>% 
    select(Type, uniprot, Thalf_wt_med, model) %>% 
    mutate(Type = if_else(Type == "TMT_light","deg","syn")) %>% 
    pivot_wider(id_cols = uniprot, names_from = c(model, Type), values_from = Thalf_wt_med) %>% 
    select(-uniprot) %>% cor(use = "pairwise.complete.obs", method = "spearman")

p_cor <- ggcorrplot::ggcorrplot(p_cor, lab=T, type = "lower")
s_cor <- ggcorrplot::ggcorrplot(s_cor, lab=T, type = "lower")

# Save data ----

# Save the data for use with the dashboard
# All of it:
dir.create(paste0(working_dir,"/script_results"))
save(list=ls(), file = paste0(working_dir,"/script_results/all_data.RData"))

# Leave just what I think I'll need:
rm(All_mods, All_mods_U,  prot_alpha_2, prot_alpha_CC_2, prot_alpha_CC_U, prot_alpha_CC_U_2, prot_alpha_U, prot_alpha_U_2, groupings, outputdir, FC_calc, filterIP2, guess_Thalf_pctHeavy, nls_mods, importIP2, match_closest, normIP2, pcnt_label_plot, proteome_mean_plot, proteome_mean_plot2, rmTMTtagIP2)

save(list=ls(), file = paste0(working_dir,"/script_results/necessary_data.RData"))
rm(list = ls())
