# Aim of this script is to run generalized linear model
# for ICD10 codes (level 2) against BP/SH/BP-SH PRS and covariates
# (Note: the analysis is restricted to non-relatives from the test sample; PRS values are standardized; only codes till Chapter XVIII selected) 

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

load("bp_prs_iid_icd10_filtered_test_nonrelatives_prs_var1.RData") # standardized PRS data
load("./icd10_iid_level_2_cbp_prev_filtered_test_nonrelatives.RData") # ICD10 data

ls()

# Start with a self-check
table(bp_f_test_nonr_var1$IID == icd_f_l2_test_nonr_prev$IID) # all TRUE

dim(bp_f_test_nonr_var1)
colnames(bp_f_test_nonr_var1)
dim(icd_f_l2_test_nonr_prev)

# Start with BP PRS; extract codes related to chapter I - XVII only
bp_test_nonr_var1_cov <- lapply(icd_f_l2_test_nonr_prev[, c(23:187)], function(x) 
		    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp, data = bp_f_test_nonr_var1, family = "binomial") })

bp_test_nonr_var1_cov <- lapply(bp_test_nonr_var1_cov, function(x) tail(summary(x)$coefficients, n = 1))

readme_bp_test_nonr_var1_cov <- "glm(binomial) results for ICD10 codes level 2, chapter I - XVII, from icd10_iid_level_2_cbp_prev_filtered_test_nonrelatives.RData against standardized BP PRS"
save(bp_test_nonr_var1_cov, readme_bp_test_nonr_var1_cov, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter_1-17_vs_bp_prs_test_nonrelatives_var1.RData") 

#rm(bp_test_nonr_var1_cov, readme_bp_test_nonr_var1_cov)

# Start with SH PRS; extract codes related to chapter I - XVII only
sh_test_nonr_var1_cov <- lapply(icd_f_l2_test_nonr_prev[, c(23:187)], function(x)
	            { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_sh, data = bp_f_test_nonr_var1, family = "binomial") })

sh_test_nonr_var1_cov <- lapply(sh_test_nonr_var1_cov, function(x) tail(summary(x)$coefficients, n = 1))

readme_sh_test_nonr_var1_cov <- "glm(binomial) results for ICD10 codes level 2, chapter I - XVII, from icd10_iid_level_2_cbp_prev_filtered_test_nonrelatives.RData against standardized SH PRS"
save(sh_test_nonr_var1_cov, readme_sh_test_nonr_var1_cov, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter_1-17_vs_sh_prs_test_nonrelatives_var1.RData")

# Start with BP-SH PRS; extract codes related to chapter I - XVII only
bp_sh_test_nonr_var1_cov <- lapply(icd_f_l2_test_nonr_prev[, c(23:187)], function(x)
	          { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp_sh, data = bp_f_test_nonr_var1, family = "binomial") })

bp_sh_test_nonr_var1_cov <- lapply(bp_sh_test_nonr_var1_cov, function(x) tail(summary(x)$coefficients, n = 1))

readme_bp_sh_test_nonr_var1_cov <- "glm(binomial) results for ICD10 codes level 2, chapter I - XVII, from icd10_iid_level_2_cbp_prev_filtered_test_nonrelatives.RData against standardized BP-SH PRS"
save(bp_sh_test_nonr_var1_cov, readme_bp_sh_test_nonr_var1_cov, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter_1-17_vs_bp_sh_prs_test_nonrelatives_var1.RData")



