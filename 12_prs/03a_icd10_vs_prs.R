# Aim of this script is to run generalized linear model
# for ICD10 codes against BP/SH/BP-SH PRS and covariates

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

load("./bp_prs_bin_iid_filtered.RData") # PRS data recoded as binary traits
load("./icd10_iid_cbp_prev_filtered.RData") # ICD10 data

ls()

# Start with a self-check
table(bp_f_bin$IID == icd_f_prev$IID) # all TRUE

dim(bp_f_bin)
colnames(bp_f_bin)
dim(icd_f_prev)

# Start with BP PRS
bp_bin_cov <- lapply(icd_f_prev[, -c(1)], function(x) 
		    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp_bin, data = bp_f_bin, family = "binomial") })

bp_bin_cov <- lapply(bp_bin_cov, function(x) tail(summary(x)$coefficients, n = 1))

readme_bp_bin_cov <- "glm() results for ICD10 codes from icd10_iid_cbp_prev_filtered.RData against binary coded BP PRS"
save(bp_bin_cov, readme_bp_bin_cov, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_vs_bp_prs.RData") 

#rm(bp_bin_cov, readme_bp_bin_cov)

# Start with SH PRS
sh_bin_cov <- lapply(icd_f_prev[, -c(1)], function(x)
	            { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_sh_bin, data = bp_f_bin, family = "binomial") })

sh_bin_cov <- lapply(sh_bin_cov, function(x) tail(summary(x)$coefficients, n = 1))

readme_sh_bin_cov <- "glm() results for ICD10 codes from icd10_iid_cbp_prev_filtered.RData against binary coded SH PRS"
save(sh_bin_cov, readme_sh_bin_cov, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_vs_sh_prs.RData")

# Start with BP-SH PRS
bp_sh_bin_cov <- lapply(icd_f_prev[, -c(1)], function(x)
	          { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp_sh_bin, data = bp_f_bin, family = "binomial") })

bp_sh_bin_cov <- lapply(bp_sh_bin_cov, function(x) tail(summary(x)$coefficients, n = 1))

readme_bp_sh_bin_cov <- "glm() results for ICD10 codes from icd10_iid_cbp_prev_filtered.RData against binary coded BP-SH PRS"
save(bp_sh_bin_cov, readme_bp_sh_bin_cov, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_vs_bp_sh_prs.RData")



