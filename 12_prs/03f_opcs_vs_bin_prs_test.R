# Aim of this script is to run generalized linear model
# for OPCS codes (not combined!!!) against SH/BP-SH PRS (binary coded; non-overlapping 01. and 0.9 quantiles of SH and BP-SH PRS) and covariates
# (Note: the analysis is restricted to non-relatives, cases from the test sample; only codes till Chapters X, Y, Z excluded) 

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

load("bp_cases_prs_iid_icd10_filtered_test_nonrelatives.RData") # load initial PRS and covariatesdata
load("binary_prs_test_nonrelatives_cases_nonoverlapping_quantiles.RData") # binary coded PRS data, non-overlapping quantiles
load("opcs_iid_cbp_prev_filtered.RData") # OPCS data, non-combined

ls()

# Reorder and self-check
ind <- match(bp_cases_f_test_nonr$IID, opcs_f_prev$IID)
opcs_f_prev <- opcs_f_prev[ind, ]
length(if_10_sh) == nrow(opcs_f_prev) # TRUE

# Select codes to exclude
chapters <- c('X', 'Y', 'Z')
excl <- lapply(chapters, function(x) grep(x, colnames(opcs_f_prev))) # test
excl <- unlist(excl)

# Start with binary coded PRS trait reflecting relatedness to 0.1 quantile of SH PRS; X, Y, Z codes excluded
if_10_sh_vs_opcs <- lapply(opcs_f_prev[, -c(1, excl)], function(x)
			                       { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_10_sh, data = bp_cases_f_test_nonr, family = "binomial") })

if_10_sh_vs_opcs <- lapply(if_10_sh_vs_opcs, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_10_sh_vs_opcs <- "glm(binomial) results for OPCS codes, X, Y, Z codes excluded, from opcs_iid_cbp_prev_filtered.RData against binary coded PRS trait reflecting relatedness to 0.1 quantile for SH"
save(if_10_sh_vs_opcs, readme_if_10_sh_vs_opcs, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_bin_prs_10_sh.RData")


# Start with binary coded PRS trait reflecting relatedness to 0.1 quantile of BP-SH PRS; X, Y, Z codes excluded
if_10_bp_sh_vs_opcs <- lapply(opcs_f_prev[, -c(1, excl)], function(x)
			                       { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_10_bp_sh, data = bp_cases_f_test_nonr, family = "binomial") })

if_10_bp_sh_vs_opcs <- lapply(if_10_bp_sh_vs_opcs, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_10_bp_sh_vs_opcs <- "glm(binomial) results for OPCS codes, X, Y, Z codes excluded, from opcs_iid_cbp_prev_filtered.RData against binary coded PRS trait reflecting relatedness to 0.1 quantile for BP-SH"
save(if_10_bp_sh_vs_opcs, readme_if_10_bp_sh_vs_opcs, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_bin_prs_10_bp_sh.RData")


# Start with binary coded PRS trait reflecting relatedness to 0.9 quantile of SH PRS; X, Y, Z codes excluded
if_90_sh_vs_opcs <- lapply(opcs_f_prev[, -c(1, excl)], function(x)
			                       { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_90_sh, data = bp_cases_f_test_nonr, family = "binomial") })

if_90_sh_vs_opcs <- lapply(if_90_sh_vs_opcs, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_90_sh_vs_opcs <- "glm(binomial) results for OPCS codes, X, Y, Z codes excluded, from opcs_iid_cbp_prev_filtered.RData against binary coded PRS trait reflecting relatedness to 0.9 quantile for SH"
save(if_90_sh_vs_opcs, readme_if_90_sh_vs_opcs, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_bin_prs_90_sh.RData")


# Start with binary coded PRS trait reflecting relatedness to 0.9 quantile of BP-SH PRS; X, Y, Z codes excluded
if_90_bp_sh_vs_opcs <- lapply(opcs_f_prev[, -c(1, excl)], function(x)
			                                                  { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_90_bp_sh, data = bp_cases_f_test_nonr, family = "binomial") })

if_90_bp_sh_vs_opcs <- lapply(if_90_bp_sh_vs_opcs, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_90_bp_sh_vs_opcs <- "glm(binomial) results for OPCS codes, X, Y, Z codes excluded, from opcs_iid_cbp_prev_filtered.RData against binary coded PRS trait reflecting relatedness to 0.9 quantile for BP-SH"
save(if_90_bp_sh_vs_opcs, readme_if_90_bp_sh_vs_opcs, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_bin_prs_90_bp_sh.RData")


