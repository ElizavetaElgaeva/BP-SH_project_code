# Aim of this script is to run generalized linear model
# for ICD10 codes (level 2) against SH/BP-SH PRS (binary coded) and covariates
# (Note: the analysis is restricted to non-relatives from the test sample; only codes till Chapter XVIII selected) 

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

load("bp_cases_prs_iid_icd10_filtered_test_nonrelatives.RData") # load initial PRS and covariatesdata
load("binary_prs_test_nonrelatives_cases.RData") # binary coded PRS data
load("./icd10_iid_cases_level_2_cbp_prev_filtered_test_nonrelatives.RData") # ICD10 data

#load("binary_prs_test_nonrelatives_cases_nonoverlapping_quantiles.RData") # binary coded PRS data; nonoverlapping quantiles
#load("./icd10_iid_cbp_prev_filtered.RData") # ICD10 data, non-combined

ls()

# Start with a self-check
length(if_10_10) == nrow(icd_cases_f_l2_test_nonr_prev) # TRUE
length(if_10_90) == nrow(icd_cases_f_l2_test_nonr_prev) # TRUE
length(if_90_10) == nrow(icd_cases_f_l2_test_nonr_prev) # TRUE
length(if_90_90) == nrow(icd_cases_f_l2_test_nonr_prev) # TRUE

#ind <- match(bp_cases_f_test_nonr$IID, icd_f_prev$IID)
#icd_f_prev <- icd_f_prev[ind, ]
#length(if_10_sh) == nrow(bp_cases_f_test_nonr) # TRUE

colnames(icd_cases_f_l2_test_nonr_prev)[1:22]
dim(icd_cases_f_l2_test_nonr_prev)

# Start with binary coded PRS trait reflecting relatedness to 0.1 quantile for both SH and BP-SH PRS; extract codes related to chapter I - XVII only
if_10_10_vs_icd <- lapply(icd_cases_f_l2_test_nonr_prev[, c(23:187)], function(x) 
		    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_10_10, data = bp_cases_f_test_nonr, family = "binomial") })

if_10_10_vs_icd <- lapply(if_10_10_vs_icd, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_10_10_vs_icd <- "glm(binomial) results for ICD10 codes level 2, chapter I - XVII, from icd10_iid_cases_level_2_cbp_prev_filtered_test_nonrelatives.RData against binary coded PRS trait reflecting relatedness to 0.1 quantile for both SH and BP-SH"
save(if_10_10_vs_icd, readme_if_10_10_vs_icd, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter_1-17_vs_bin_prs_10_10.RData") 


#if_10_sh_vs_icd <- lapply(icd_f_prev[, -c(1)], function(x)
#			                      { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_10_sh, data = bp_cases_f_test_nonr, family = "binomial") })

#if_10_sh_vs_icd <- lapply(if_10_sh_vs_icd, function(x) tail(summary(x)$coefficients, n = 1))

#readme_if_10_sh_vs_icd <- "glm(binomial) results for ICD10 codes, from icd10_iid_cbp_prev_filtered.RData against binary coded PRS trait reflecting relatedness to 0.1 quantile for SH"
#save(if_10_sh_vs_icd, readme_if_10_sh_vs_icd, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_vs_bin_prs_10_sh.RData")


# Start with binary coded PRS trait reflecting relatedness to 0.1 quantile for SH and 0.9 for BP-SH PRS; extract codes related to chapter I - XVII only
if_10_90_vs_icd <- lapply(icd_cases_f_l2_test_nonr_prev[, c(23:187)], function(x)
	            { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_10_90, data = bp_cases_f_test_nonr, family = "binomial") })

if_10_90_vs_icd <- lapply(if_10_90_vs_icd, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_10_90_vs_icd <- "glm(binomial) results for ICD10 codes level 2, chapter I - XVII, from icd10_iid_cases_level_2_cbp_prev_filtered_test_nonrelatives.RData against binary coded PRS trait reflecting relatedness to 0.1 quantile for SH and 0.9 for BP-SH PRS"
save(if_10_90_vs_icd, readme_if_10_90_vs_icd, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter_1-17_vs_bin_prs_10_90.RData")


#if_10_bp_sh_vs_icd <- lapply(icd_f_prev[, -c(1)], function(x)
#			                      { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_10_bp_sh, data = bp_cases_f_test_nonr, family = "binomial") })

#if_10_bp_sh_vs_icd <- lapply(if_10_bp_sh_vs_icd, function(x) tail(summary(x)$coefficients, n = 1))

#readme_if_10_bp_sh_vs_icd <- "glm(binomial) results for ICD10 codes, from icd10_iid_cbp_prev_filtered.RData against binary coded PRS trait reflecting relatedness to 0.1 quantile for BP-SH PRS"
#save(if_10_bp_sh_vs_icd, readme_if_10_bp_sh_vs_icd, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_vs_bin_prs_10_bp_sh.RData")


# Start with binary coded PRS trait reflecting relatedness to 0.9 quantile for SH and 0.1 for BP-SH PRS; extract codes related to chapter I - XVII only
if_90_10_vs_icd <- lapply(icd_cases_f_l2_test_nonr_prev[, c(23:187)], function(x)
	          { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_90_10, data = bp_cases_f_test_nonr, family = "binomial") })

if_90_10_vs_icd <- lapply(if_90_10_vs_icd, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_90_10_vs_icd <- "glm(binomial) results for ICD10 codes level 2, chapter I - XVII, from icd10_iid_cases_level_2_cbp_prev_filtered_test_nonrelatives.RData against binary coded PRS trait reflecting relatedness to 0.9 quantile for SH and 0.1 for BP-SH PRS"
save(if_90_10_vs_icd, readme_if_90_10_vs_icd, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter_1-17_vs_bin_sh_prs_90_10.RData")


#if_90_sh_vs_icd <- lapply(icd_f_prev[, -c(1)], function(x)
#			                    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_90_sh, data = bp_cases_f_test_nonr, family = "binomial") })

#if_90_sh_vs_icd <- lapply(if_90_sh_vs_icd, function(x) tail(summary(x)$coefficients, n = 1))

#readme_if_90_sh_vs_icd <- "glm(binomial) results for ICD10 codes, from icd10_iid_cbp_prev_filtered.RData against binary coded PRS trait reflecting relatedness to 0.9 quantile for SH PRS"
#save(if_90_sh_vs_icd, readme_if_90_sh_vs_icd, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_vs_bin_sh_prs_90_sh.RData")


# Start with binary coded PRS trait reflecting relatedness to 0.9 quantile for both SH and BP-SH PRS; extract codes related to chapter I - XVII only
if_90_90_vs_icd <- lapply(icd_cases_f_l2_test_nonr_prev[, c(23:187)], function(x)
			                      { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_90_90, data = bp_cases_f_test_nonr, family = "binomial") })

if_90_90_vs_icd <- lapply(if_90_90_vs_icd, function(x) tail(summary(x)$coefficients, n = 1))

readme_if_90_90_vs_icd <- "glm(binomial) results for ICD10 codes level 2, chapter I - XVII, from icd10_iid_cases_level_2_cbp_prev_filtered_test_nonrelatives.RData against binary coded PRS trait reflecting relatedness to 0.9 quantile for both SH and BP-SH"
save(if_90_90_vs_icd, readme_if_90_90_vs_icd, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter_1-17_vs_bin_prs_90_90.RData")


#if_90_bp_sh_vs_icd <- lapply(icd_f_prev[, -c(1)], function(x)
#			                                                { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + if_90_bp_sh, data = bp_cases_f_test_nonr, family = "binomial") })

#if_90_bp_sh_vs_icd <- lapply(if_90_bp_sh_vs_icd, function(x) tail(summary(x)$coefficients, n = 1))

#readme_if_90_bp_sh_vs_icd <- "glm(binomial) results for ICD10 codes, from icd10_iid_cbp_prev_filtered.RData against binary coded PRS trait reflecting relatedness to 0.9 quantile for BP-SH"
#save(if_90_bp_sh_vs_icd, readme_if_90_bp_sh_vs_icd, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_vs_bin_prs_90_bp_sh.RData")



