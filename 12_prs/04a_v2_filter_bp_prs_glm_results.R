# Aim of this script is to filter glm() results for ICD10 and OPCS codes
# against standardized BP PRS (test sample non-relatieves only)

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Start with BP PRS
load("./icd10_vs_bp_prs_test_nonrelatives_var1.RData")
load("./opcs_vs_bp_prs_test_nonrelatives_var1.RData")
ls()

thr <- 0.05/(3*(length(bp_test_nonr_var1_cov) + length(opcs_bp_test_nonr_var1))) # p = 3.646973e-05 is a significance threshold for glm() analyses

# ICD10
pval_icd <- lapply(bp_test_nonr_var1_cov, function(x) x[, 4])
pval_icd_i <- which(pval_icd < thr)

or_icd <- lapply(bp_test_nonr_var1_cov, function(x) exp(x[, 1]))
summary(unlist(or_icd))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9229  1.0606  1.0953  1.1013  1.1382  1.3870
beta_icd <- lapply(bp_test_nonr_var1_cov, function(x) abs(x[, 1]))
or_icd_i <- which(beta_icd > log(2))

bp_icd_i <- intersect(pval_icd_i, or_icd_i) # 0 elemets
#bp_icd <- bp_test_nonr_var1_cov[bp_icd_i] # glm() results for ICD10 codes against BP PRS, significant and abs(beta) > log(2)
#or_bp_icd <- or_icd[bp_icd_i]
bp_icd <- bp_test_nonr_var1_cov[pval_icd_i] # glm() results for ICD10 codes against BP PRS, significant
or_bp_icd <- or_icd[pval_icd_i]

readme_bp_icd <- "glm() results for ICD10 codes against BP PRS, significant"
save(readme_bp_icd, or_bp_icd, bp_icd, file = "./glm_icd_bp_prs_filtered_test_nonrelatives.RData")

# OPCS
pval_opcs <- lapply(opcs_bp_test_nonr_var1, function(x) x[, 4])
pval_opcs_i <- which(pval_opcs < thr)

or_opcs <- lapply(opcs_bp_test_nonr_var1, function(x) exp(x[, 1]))
summary(unlist(or_opcs))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9305  1.0440  1.0818  1.0847  1.1164  1.4089
beta_opcs <- lapply(opcs_bp_test_nonr_var1, function(x) abs(x[, 1]))
or_opcs_i <- which(beta_opcs > log(2))

bp_opcs_i <- intersect(pval_opcs_i, or_opcs_i) # 0 elements
#bp_opcs <- opcs_bp_test_nonr_var1[bp_opcs_i] # glm() results for ICD10 codes against BP PRS, significant and abs(beta) > log(2)
#or_bp_opcs <- or_opcs[bp_opcs_i]

bp_opcs <- opcs_bp_test_nonr_var1[pval_opcs_i] # glm() results for ICD10 codes against BP PRS, significant
or_bp_opcs <- or_opcs[pval_opcs_i]

readme_bp_opcs <- "glm() results for OPCS codes against BP PRS, significant"
save(readme_bp_opcs, or_bp_opcs, bp_opcs, file = "./glm_opcs_bp_prs_filtered_test_nonrelatives.RData")



