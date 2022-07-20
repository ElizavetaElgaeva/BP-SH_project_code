# Aim of this script is to run generalized linear model
# for ICD10 codes (not combined!!!) against SH/BP-SH PRS (binary coded; non-overlapping 0.1 and 0.9 quantiles of BP-SH and SH PRS) and covariates
# (Note: the analysis is restricted to non-relatives, cases from the test sample; only codes till Chapter XVIII selected) 

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

load("bp_cases_prs_iid_icd10_filtered_test_nonrelatives.RData") # load initial PRS and covariates data; test sample
load("binary_prs_test_nonrelatives_cases_nonoverlapping_quantiles.RData") # binary coded PRS data; non-overlapping 0.1 and 0.9 BP-SH and SH PRS quantiles
load("./icd10_iid_cbp_prev_filtered.RData") # ICD10 data, non-combined

ls()

# Reorder abd start with a self-check
ind <- match(bp_cases_f_test_nonr$IID, icd_f_prev$IID)
icd_f_prev <- icd_f_prev[ind, ]
length(if_10_sh) == nrow(bp_cases_f_test_nonr) # TRUE

# Select codes to exclude
chapters <- c('R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', 'U')
excl <- lapply(chapters, function(x) grep(x, colnames(icd_f_prev))) # test
excl <- unlist(excl)

# Start with binary coded PRS trait reflecting relatedness to 0.1 quantile of SH PRS; extract codes related to chapter I - XVII only
if_10_sh_vs_icd <- lapply(icd_f_prev[, -c(1, excl)], function(x)
			                      { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_10_sh, data = bp_cases_f_test_nonr, family = "binomial") })

if_10_sh_vs_icd <- lapply(if_10_sh_vs_icd, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_10_sh_vs_icd <- "glm(binomial) results for ICD10 codes, chapter I - XVII, from icd10_iid_cbp_prev_filtered.RData against binary coded PRS trait reflecting relatedness to 0.1 quantile for SH"
save(if_10_sh_vs_icd, readme_if_10_sh_vs_icd, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_chapter1-17_vs_bin_prs_10_sh.RData")


# Start with binary coded PRS trait reflecting relatedness to 0.1 quantile of BP-SH PRS; extract codes related to chapter I - XVII only
if_10_bp_sh_vs_icd <- lapply(icd_f_prev[, -c(1, excl)], function(x)
			                      { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_10_bp_sh, data = bp_cases_f_test_nonr, family = "binomial") })

if_10_bp_sh_vs_icd <- lapply(if_10_bp_sh_vs_icd, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_10_bp_sh_vs_icd <- "glm(binomial) results for ICD10 codes, chapter I - XVII, from icd10_iid_cbp_prev_filtered.RData against binary coded PRS trait reflecting relatedness to 0.1 quantile for BP-SH PRS"
save(if_10_bp_sh_vs_icd, readme_if_10_bp_sh_vs_icd, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_chapter1-17_vs_bin_prs_10_bp_sh.RData")


# Start with binary coded PRS trait reflecting relatedness to 0.9 quantile of SH PRS; extract codes related to chapter I - XVII only
if_90_sh_vs_icd <- lapply(icd_f_prev[, -c(1, excl)], function(x)
			                    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_90_sh, data = bp_cases_f_test_nonr, family = "binomial") })

if_90_sh_vs_icd <- lapply(if_90_sh_vs_icd, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_90_sh_vs_icd <- "glm(binomial) results for ICD10 codes, chapter I - XVII, from icd10_iid_cbp_prev_filtered.RData against binary coded PRS trait reflecting relatedness to 0.9 quantile for SH PRS"
save(if_90_sh_vs_icd, readme_if_90_sh_vs_icd, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_chapter1-17vs_bin_sh_prs_90_sh.RData")


# Start with binary coded PRS trait reflecting relatedness to 0.9 quantile of BP-SH PRS; extract codes related to chapter I - XVII only
if_90_bp_sh_vs_icd <- lapply(icd_f_prev[, -c(1, excl)], function(x)
			                                                { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_90_bp_sh, data = bp_cases_f_test_nonr, family = "binomial") })

if_90_bp_sh_vs_icd <- lapply(if_90_bp_sh_vs_icd, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_90_bp_sh_vs_icd <- "glm(binomial) results for ICD10 codes, chapter I - XVII, from icd10_iid_cbp_prev_filtered.RData against binary coded PRS trait reflecting relatedness to 0.9 quantile for BP-SH"
save(if_90_bp_sh_vs_icd, readme_if_90_bp_sh_vs_icd, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_chapter1-17_vs_bin_prs_90_bp_sh.RData")



