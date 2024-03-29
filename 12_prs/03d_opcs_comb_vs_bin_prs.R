# Aim of this script is to run generalized linear model
# for OPCS codes (combined to the level 2!!!) against SH/BP-SH PRS (binary coded; overlapping 0.1 / 0.9 quantiles of SH and BP-SH PRS) and covariates
# (Note: the analysis is restricted to non-relatives, cases from the test sample; only codes till Chapters X, Y, Z excluded) 

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

load("bp_cases_prs_iid_icd10_filtered_test_nonrelatives.RData") # load initial PRS and covariatesdata; test sample
load("binary_prs_test_nonrelatives_cases.RData") # binary coded PRS data; overlapping quantiles
load("./opcs_iid_cases_level_2_cbp_prev_filtered_test_nonrelatives.RData") # OPCS data; combined

ls()

# Start with a self-check
length(if_10_10) == nrow(opcs_cases_f_l2_test_nonr_prev) # TRUE
length(if_10_90) == nrow(opcs_cases_f_l2_test_nonr_prev) # TRUE
length(if_90_10) == nrow(opcs_cases_f_l2_test_nonr_prev) # TRUE
length(if_90_90) == nrow(opcs_cases_f_l2_test_nonr_prev) # TRUE

colnames(opcs_cases_f_l2_test_nonr_prev)[1:22]
dim(opcs_cases_f_l2_test_nonr_prev)

# Start with binary coded PRS trait reflecting relatedness to 0.1 quantile for both SH and BP-SH PRS; X, Y, Z codes excluded
if_10_10_vs_opcs <- lapply(opcs_cases_f_l2_test_nonr_prev[, c(23:154)], function(x) 
		    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_10_10, data = bp_cases_f_test_nonr, family = "binomial") })

if_10_10_vs_opcs <- lapply(if_10_10_vs_opcs, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_10_10_vs_opcs <- "glm(binomial) results for OPCS codes level 2, X, Y, Z codes excluded, from opcs_iid_cases_level_2_cbp_prev_filtered_test_nonrelatives.RData against binary coded PRS trait reflecting relatedness to 0.1 quantile for both SH and BP-SH"
save(if_10_10_vs_opcs, readme_if_10_10_vs_opcs, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_level_2_no_xyz_vs_bin_prs_10_10.RData") 

# Start with binary coded PRS trait reflecting relatedness to 0.1 quantile for SH and 0.9 for BP-SH PRS; X, Y, Z codes excluded
if_10_90_vs_opcs <- lapply(opcs_cases_f_l2_test_nonr_prev[, c(23:154)], function(x)
	            { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_10_90, data = bp_cases_f_test_nonr, family = "binomial") })

if_10_90_vs_opcs <- lapply(if_10_90_vs_opcs, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_10_90_vs_opcs <- "glm(binomial) results for ICD10 codes level 2,  X, Y, Z codes excluded, from opcs_iid_cases_level_2_cbp_prev_filtered_test_nonrelatives.RData against binary coded PRS trait reflecting relatedness to 0.1 quantile for SH and 0.9 for BP-SH PRS"
save(if_10_90_vs_opcs, readme_if_10_90_vs_opcs, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_level_2_no_xyz_vs_bin_prs_10_90.RData")

# Start with binary coded PRS trait reflecting relatedness to 0.9 quantile for SH and 0.1 for BP-SH PRS; X, Y, Z codes excluded
if_90_10_vs_opcs <- lapply(opcs_cases_f_l2_test_nonr_prev[, c(23:154)], function(x)
	          { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_90_10, data = bp_cases_f_test_nonr, family = "binomial") })

if_90_10_vs_opcs <- lapply(if_90_10_vs_opcs, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_90_10_vs_opcs <- "glm(binomial) results for ICD10 codes level 2, X, Y, Z codes excluded, from opcs_iid_cases_level_2_cbp_prev_filtered_test_nonrelatives.RData against binary coded PRS trait reflecting relatedness to 0.9 quantile for SH and 0.1 for BP-SH PRS"
save(if_90_10_vs_opcs, readme_if_90_10_vs_opcs, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_level_2_no_xyz_vs_bin_sh_prs_90_10.RData")

# Start with binary coded PRS trait reflecting relatedness to 0.9 quantile for both SH and BP-SH PRS; X, Y, Z codes excluded
if_90_90_vs_opcs <- lapply(opcs_cases_f_l2_test_nonr_prev[, c(23:154)], function(x)
			                      { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_90_90, data = bp_cases_f_test_nonr, family = "binomial") })

if_90_90_vs_opcs <- lapply(if_90_90_vs_opcs, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_90_90_vs_opcs <- "glm(binomial) results for ICD10 codes level 2, X, Y, Z codes excluded, from opcs_iid_cases_level_2_cbp_prev_filtered_test_nonrelatives.RData against binary coded PRS trait reflecting relatedness to 0.9 quantile for both SH and BP-SH"
save(if_90_90_vs_opcs, readme_if_90_90_vs_opcs, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_level_2_no_xyz_vs_bin_prs_90_90.RData")




