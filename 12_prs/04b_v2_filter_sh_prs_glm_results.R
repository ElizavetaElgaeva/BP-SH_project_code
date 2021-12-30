# Aim of this script is to filter glm() results for ICD10 and OPCS codes
# against standardized SH PRS (test sample non-relatieves only)

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Start with SH PRS
load("./icd10_vs_sh_prs_test_nonrelatives_var1.RData")
load("./opcs_vs_sh_prs_test_nonrelatives_var1.RData")
ls()

thr <- 0.05/(3*(length(sh_test_nonr_var1_cov) + length(opcs_sh_test_nonr_var1))) # p = 3.646973e-05 is a significance threshold for glm() analyses

# ICD10
pval_icd <- lapply(sh_test_nonr_var1_cov, function(x) x[, 4])
pval_icd_i <- which(pval_icd < thr)

or_icd <- lapply(sh_test_nonr_var1_cov, function(x) exp(x[, 1]))
summary(unlist(or_icd))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.8985  1.0892  1.1322  1.1410  1.1878  1.4220
beta_icd <- lapply(sh_test_nonr_var1_cov, function(x) abs(x[, 1]))
or_icd_i <- which(beta_icd > log(2))

sh_icd_i <- intersect(pval_icd_i, or_icd_i) # 0 elemets
#sh_icd <- bp_test_nonr_var1_cov[bp_icd_i] # glm() results for ICD10 codes against SH PRS, significant and abs(beta) > log(2)
#or_sh_icd <- or_icd[sh_icd_i]
sh_icd <- sh_test_nonr_var1_cov[pval_icd_i] # glm() results for ICD10 codes against SH PRS, significant
or_sh_icd <- or_icd[pval_icd_i]

readme_sh_icd <- "glm() results for ICD10 codes against SH PRS, significant"
save(readme_sh_icd, or_sh_icd, sh_icd, file = "./glm_icd_sh_prs_filtered_test_nonrelatives.RData")

# OPCS
pval_opcs <- lapply(opcs_sh_test_nonr_var1, function(x) x[, 4])
pval_opcs_i <- which(pval_opcs < thr)

or_opcs <- lapply(opcs_sh_test_nonr_var1, function(x) exp(x[, 1]))
summary(unlist(or_opcs))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.903   1.057   1.107   1.113   1.163   1.396
beta_opcs <- lapply(opcs_sh_test_nonr_var1, function(x) abs(x[, 1]))
or_opcs_i <- which(beta_opcs > log(2))

sh_opcs_i <- intersect(pval_opcs_i, or_opcs_i) # 0 elements
#sh_opcs <- opcs_sh_test_nonr_var1[sh_opcs_i] # glm() results for ICD10 codes against SH PRS, significant and abs(beta) > log(2)
#or_sh_opcs <- or_opcs[sh_opcs_i]

sh_opcs <- opcs_sh_test_nonr_var1[pval_opcs_i] # glm() results for ICD10 codes against SH PRS, significant
or_sh_opcs <- or_opcs[pval_opcs_i]

readme_sh_opcs <- "glm() results for OPCS codes against SH PRS, significant"
save(readme_sh_opcs, or_sh_opcs, sh_opcs, file = "./glm_opcs_sh_prs_filtered_test_nonrelatives.RData")



