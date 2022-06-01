# Aim of this script is to run generalized linear model
# for OPCS codes (level 2) against BP/SH/BP-SH PRS and covariates
# (Note: the analysis is restricted to non-relatives from the test sample; PRS values are standardized; X, Y, Z codes excluded)

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

load("bp_prs_iid_icd10_filtered_test_nonrelatives_prs_var1.RData") # standardized PRS data
load("./opcs_iid_level_2_cbp_prev_filtered_test_nonrelatives.RData") # OPCS data

ls()

# Start with a self-check
table(bp_f_test_nonr_var1$IID == opcs_f_l2_test_nonr_prev$IID) # all TRUE

dim(bp_f_test_nonr_var1)
colnames(bp_f_test_nonr_var1)
dim(opcs_f_l2_test_nonr_prev)

# Start with BP PRS; X, Y, Z codes excluded
opcs_bp_test_nonr_var1 <- lapply(opcs_f_l2_test_nonr_prev[, c(23:154)], function(x) 
		    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp, data = bp_f_test_nonr_var1, family = "binomial") })

opcs_bp_test_nonr_var1 <- lapply(opcs_bp_test_nonr_var1, function(x) tail(summary(x)$coefficients, n = 1))

readme_opcs_bp_test_nonr_var1 <- "glm(binomial) results for OPCS codes (level 2; X, Y, Z codes excluded) from opcs_iid_level_2_cbp_prev_filtered_test_nonrelatives.RData against standardized BP PRS"
save(opcs_bp_test_nonr_var1, readme_opcs_bp_test_nonr_var1, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_level_2_no_xyz_vs_bp_prs_test_nonrelatives_var1.RData") 

# Start with SH PRS; X, Y, Z codes excluded
opcs_sh_test_nonr_var1 <- lapply(opcs_f_l2_test_nonr_prev[, c(23:154)], function(x)
	            { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_sh, data = bp_f_test_nonr_var1, family = "binomial") })

opcs_sh_test_nonr_var1 <- lapply(opcs_sh_test_nonr_var1, function(x) tail(summary(x)$coefficients, n = 1))

readme_opcs_sh_test_nonr_var1 <- "glm(binomial) results for OPCS codes (level 2; X, Y, Z codes excluded) from opcs_iid_level_2_cbp_prev_filtered_test_nonrelatives.RData against standardized SH PRS"
save(opcs_sh_test_nonr_var1, readme_opcs_sh_test_nonr_var1, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_level_2_no_xyz_vs_sh_prs_test_nonrelatives_var1.RData")

# Start with BP-SH PRS; X, Y, Z codes excluded
opcs_bp_sh_test_nonr_var1 <- lapply(opcs_f_l2_test_nonr_prev[, c(23:154)], function(x)
	          { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp_sh, data = bp_f_test_nonr_var1, family = "binomial") })

opcs_bp_sh_test_nonr_var1 <- lapply(opcs_bp_sh_test_nonr_var1, function(x) tail(summary(x)$coefficients, n = 1))

readme_opcs_bp_sh_test_nonr_var1 <- "glm(binomial) results for OPCS codes (level 2; X, Y, Z codes excluded) from opcs_iid_level_2_cbp_prev_filtered_test_nonrelatives.RData against standardized BP-SH PRS"
save(opcs_bp_sh_test_nonr_var1, readme_opcs_bp_sh_test_nonr_var1, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_level_2_no_xyz_vs_bp_sh_prs_test_nonrelatives_var1.RData")



