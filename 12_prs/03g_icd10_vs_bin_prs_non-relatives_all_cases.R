# Aim of this script is to run generalized linear model
# for ICD10 codes (combined to the level 2) against SH/BP-SH PRS (binary coded; non-overlapping 0.1 and 0.9 quantiles of SH and BP-SH PRS) and covariates
# (Note: the analysis is restricted to non-relatives, cases from the whole sample; only codes till Chapter XVIII selected) 

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

load("bp_cases_prs_iid_icd10_filtered_all_nonrel.RData") # load initial PRS and covariatesdata; all non-relatives cases
load("binary_prs_all_nonrel_cases_nonoverlapping_quantiles.RData") # binary coded PRS data; nonoverlapping quantiles; all non-relatives cases
load("icd10_level_2_iid_cases_cbp_prev_filtered_all_nonrel.RData") # ICD10 data combined to the level 2; all non-relatives cases

ls()

# Start with a self-check
length(if_10_sh) == nrow(icd_l2_cases_nonr_prev) # TRUE
length(if_10_bp_sh) == nrow(icd_l2_cases_nonr_prev) # TRUE
length(if_90_bp_sh) == nrow(icd_l2_cases_nonr_prev) # TRUE
length(if_90_sh) == nrow(icd_l2_cases_nonr_prev) # TRUE

# Select codes to exclude
chapters <- c('R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', 'U')
excl <- lapply(chapters, function(x) grep(x, colnames(icd_l2_cases_nonr_prev)))
excl <- unlist(excl)
excl <- unique(c(1:22, excl)) # exclude firts 22 columns since they do not contain codes at all

# Start with binary coded PRS trait reflecting relatedness to 0.1 quantile of SH PRS; extract codes related to chapter I - XVII only
if_10_sh_vs_icd_l2 <- lapply(icd_l2_cases_nonr_prev[, -excl], function(x)
                                             { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_10_sh, data = bp_cases_nonr, family = "binomial") })

if_10_sh_vs_icd_l2 <- lapply(if_10_sh_vs_icd_l2, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_10_sh_vs_icd_l2 <- "glm(binomial) results for ICD10 codes, level 2, chapter I - XVII, from icd10_level_2_iid_cases_cbp_prev_filtered_all_nonrel.RData against binary coded PRS trait reflecting relatedness to 0.1 quantile for SH"
save(if_10_sh_vs_icd_l2, readme_if_10_sh_vs_icd_l2, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter1-17_vs_bin_prs_10_sh_all_cases.RData")

# Start with binary coded PRS trait reflecting relatedness to 0.1 quantile of BP-SH PRS; extract codes related to chapter I - XVII only
if_10_bp_sh_vs_icd_l2 <- lapply(icd_l2_cases_nonr_prev[, -excl], function(x)
                                             { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_10_bp_sh, data = bp_cases_nonr, family = "binomial") })

if_10_bp_sh_vs_icd_l2 <- lapply(if_10_bp_sh_vs_icd_l2, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_10_bp_sh_vs_icd_l2 <- "glm(binomial) results for ICD10 codes, level 2, chapter I - XVII, from icd10_level_2_iid_cases_cbp_prev_filtered_all_nonrel.RData against binary coded PRS trait reflecting relatedness to 0.1 quantile for BP-SH PRS"
save(if_10_bp_sh_vs_icd_l2, readme_if_10_bp_sh_vs_icd_l2, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter1-17_vs_bin_prs_10_bp_sh_all_cases.RData")

# Start with binary coded PRS trait reflecting relatedness to 0.9 quantile of SH PRS; extract codes related to chapter I - XVII only
if_90_sh_vs_icd_l2 <- lapply(icd_l2_cases_nonr_prev[, -excl], function(x)
                                           { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_90_sh, data = bp_cases_nonr, family = "binomial") })

if_90_sh_vs_icd_l2 <- lapply(if_90_sh_vs_icd_l2, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_90_sh_vs_icd_l2 <- "glm(binomial) results for ICD10 codes, level 2, chapter I - XVII, from icd10_level_2_iid_cases_cbp_prev_filtered_all_nonrel.RData against binary coded PRS trait reflecting relatedness to 0.9 quantile for SH PRS"
save(if_90_sh_vs_icd_l2, readme_if_90_sh_vs_icd_l2, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter1-17_vs_bin_sh_prs_90_sh_all_cases.RData")


# Start with binary coded PRS trait reflecting relatedness to 0.9 quantile of BP-SH PRS; extract codes related to chapter I - XVII only
if_90_bp_sh_vs_icd_l2 <- lapply(icd_l2_cases_nonr_prev[, -excl], function(x)
                                                                       { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_90_bp_sh, data = bp_cases_nonr, family = "binomial") })

if_90_bp_sh_vs_icd_l2 <- lapply(if_90_bp_sh_vs_icd_l2, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_90_bp_sh_vs_icd_l2 <- "glm(binomial) results for ICD10 codes, level 2, chapter I - XVII, from icd10_level_2_iid_cases_cbp_prev_filtered_all_nonrel.RData against binary coded PRS trait reflecting relatedness to 0.9 quantile for BP-SH"
save(if_90_bp_sh_vs_icd_l2, readme_if_90_bp_sh_vs_icd_l2, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter1-17_vs_bin_prs_90_bp_sh_all_cases.RData")


