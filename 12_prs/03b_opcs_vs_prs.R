# Aim of this script is to run generalized linear model
# for OPCS codes against BP/SH/BP-SH PRS and covariates

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

load("./bp_prs_bin_iid_filtered.RData") # PRS data recoded as binary traits
load("./opcs_iid_cbp_prev_filtered.RData") # OPCS data

ls()

# Start with a self-check
table(bp_f_bin$IID == opcs_f_prev$IID) # all TRUE

dim(bp_f_bin)
colnames(bp_f_bin)
dim(opcs_f_prev)

# Start with BP PRS
opcs_bp_bin <- lapply(opcs_f_prev[, -c(1)], function(x) 
		    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp_bin, data = bp_f_bin, family = "binomial") })

opcs_bp_bin <- lapply(opcs_bp_bin, function(x) tail(summary(x)$coefficients, n = 1))

readme_opcs_bp_bin <- "glm() results for OPCS codes from opcs_iid_cbp_prev_filtered.RData against binary coded BP PRS"
save(opcs_bp_bin, readme_opcs_bp_bin, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_vs_bp_prs.RData") 

# Start with SH PRS
opcs_sh_bin <- lapply(opcs_f_prev[, -c(1)], function(x)
	            { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_sh_bin, data = bp_f_bin, family = "binomial") })

opcs_sh_bin <- lapply(opcs_sh_bin, function(x) tail(summary(x)$coefficients, n = 1))

readme_opcs_sh_bin <- "glm() results for OPCS codes from opcs_iid_cbp_prev_filtered.RData against binary coded SH PRS"
save(opcs_sh_bin, readme_opcs_sh_bin, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_vs_sh_prs.RData")

# Start with BP-SH PRS
opcs_bp_sh_bin <- lapply(opcs_f_prev[, -c(1)], function(x)
	          { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp_sh_bin, data = bp_f_bin, family = "binomial") })

opcs_bp_sh_bin <- lapply(opcs_bp_sh_bin, function(x) tail(summary(x)$coefficients, n = 1))

readme_opcs_bp_sh_bin <- "glm() results for OPCS codes from opcs_iid_cbp_prev_filtered.RData against binary coded BP-SH PRS"
save(opcs_bp_sh_bin, readme_opcs_bp_sh_bin, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_vs_bp_sh_prs.RData")



