# Aim of this script is to standardize PRS values for CBP UGIT and SGIT in non-relatives cases (whole sample) and recode them to binary traits

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

# Load PRS data (whole sample, non-relatieves, cases only)
load("bp_cases_prs_iid_icd10_filtered_all_nonrel.RData")
ls()

# Standardize PRS
colnames(bp_cases_nonr)
sh <- bp_cases_nonr$"prs_sh"/sqrt(var(bp_cases_nonr$"prs_sh"))
bp_sh <- bp_cases_nonr$"prs_bp_sh"/sqrt(var(bp_cases_nonr$"prs_bp_sh"))

# Recode PRS
l_sh <- rep(0, length(sh)) # new variable derived from sh
l_sh[sh >= quantile(sh, 0.9)] <- 90 # recode values from 0.9 quantile 
l_sh[sh <= quantile(sh, 0.1)] <- 10 # recode values from 0.1 quantile
table(l_sh)

l_bp_sh <- rep(0, length(bp_sh)) # new variable derived from bp_sh
l_bp_sh[bp_sh >= quantile(bp_sh, 0.9)] <- 90 # recode values from 0.9 quantile
l_bp_sh[bp_sh <= quantile(bp_sh, 0.1)] <- 10 # recode values from 0.1 quantile
table(l_bp_sh)

table(l_sh,l_bp_sh)

# Recode PRS as binary trait based on relatedness to 0.1/0.9 quantile for CBP UGIT (BP-SH) and SGIT (SH) 
if_10_sh <- if_10_bp_sh <- if_90_sh <- if_90_bp_sh <- rep(0, length(sh)) # generate new blank variables reflecting relatedness to specific group

if_10_sh[which(l_sh == 10)] <- 1 # in 0.1 quantile for SH
table(if_10_sh)

if_10_bp_sh[which(l_bp_sh == 10)] <- 1 # in 0.1 quantile for BP-SH
table(if_10_bp_sh)

if_90_sh[which(l_sh == 90)] <- 1 # in 0.9 quantile for SH
table(if_90_sh)

if_90_bp_sh[which(l_bp_sh == 90)] <- 1 # in 0.9 quantile for BP-SH
table(if_90_bp_sh)

save(if_10_sh, if_10_bp_sh, if_90_sh, if_90_bp_sh, file = "binary_prs_all_nonrel_cases_nonoverlapping_quantiles.RData")


