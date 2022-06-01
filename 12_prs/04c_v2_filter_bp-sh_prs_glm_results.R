# Aim of this script is to filter glm() results for ICD10 and OPCS codes (both combined to level 2)
# against standardized BP-SH PRS (test sample non-relatieves only)

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Start with BP-SH PRS
load("./icd10_level_2_vs_bp_sh_prs_test_nonrelatives_var1.RData")
load("./opcs_level_2_vs_bp_sh_prs_test_nonrelatives_var1.RData")
ls()

thr <- 0.05/(3*(length(bp_sh_test_nonr_var1_cov) + length(opcs_bp_sh_test_nonr_var1) + 6)) # p = 3.687316e-05 is a significance threshold for glm() analyses; 6 - is a number of pain traits, glm results for them were included later

# ICD10
pval_icd <- lapply(bp_sh_test_nonr_var1_cov, function(x) x[, 4])
pval_icd_i <- which(pval_icd < thr) # 0 elements

or_icd <- lapply(bp_sh_test_nonr_var1_cov, function(x) exp(x[, 1]))
summary(unlist(or_icd))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.8955  0.9802  0.9929  0.9953  1.0111  1.1093
beta_icd <- lapply(bp_sh_test_nonr_var1_cov, function(x) abs(x[, 1]))
or_icd_i <- which(beta_icd > log(2))

bp_sh_icd_i <- intersect(pval_icd_i, or_icd_i) # 0 elemets
#bp_sh_icd <- bp_sh_test_nonr_var1_cov[bp_sh_icd_i] # glm() results for ICD10 codes against BP-SH PRS, significant and abs(beta) > log(2)
#or_bp_sh_icd <- or_icd[bp_sh_icd_i]
#bp_sh_icd <- bp_sh_test_nonr_var1_cov[pval_icd_i] # glm() results for ICD10 codes against BP-SH PRS, significant; N = 0
#or_bp_sh_icd <- or_icd[pval_icd_i]

#readme_bp_sh_icd <- "glm() results for ICD10 codes level 2 against BP PRS, significant"
#save(readme_bp_sh_icd, or_bp_sh_icd, bp_sh_icd, file = "./glm_icd_level_2_bp_sh_prs_filtered_test_nonrelatives.RData")

# OPCS
pval_opcs <- lapply(opcs_bp_sh_test_nonr_var1, function(x) x[, 4])
pval_opcs_i <- which(pval_opcs < thr)

or_opcs <- lapply(opcs_bp_sh_test_nonr_var1, function(x) exp(x[, 1]))
summary(unlist(or_opcs))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.8920  0.9782  0.9996  0.9998  1.0201  1.1563
beta_opcs <- lapply(opcs_bp_sh_test_nonr_var1, function(x) abs(x[, 1]))
or_opcs_i <- which(beta_opcs > log(2))

bp_sh_opcs_i <- intersect(pval_opcs_i, or_opcs_i) # 0 elements
#bp_sh_opcs <- opcs_bp_sh_test_nonr_var1[bp_sh_opcs_i] # glm() results for ICD10 codes against BP-SH PRS, significant and abs(beta) > log(2)
#or_bp_sh_opcs <- or_opcs[bp_sh_opcs_i]

bp_sh_opcs <- opcs_bp_sh_test_nonr_var1[pval_opcs_i] # glm() results for ICD10 codes against BP-SH PRS, significant; N = 1
or_bp_sh_opcs <- or_opcs[pval_opcs_i]

readme_bp_sh_opcs <- "glm() results for OPCS codes level 2 against BP-SH PRS, significant"
save(readme_bp_sh_opcs, or_bp_sh_opcs, bp_sh_opcs, file = "./glm_opcs_level_2_bp_sh_prs_filtered_test_nonrelatives.RData")



