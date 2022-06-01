# Aim of this script is to filter glm() results for ICD10 and OPCS codes (both combined to level 2)
# against standardized SH PRS (test sample non-relatieves only)
# Note: ICD10 from chapter I - XVII only; X, Y, Z OPCS codes excluded

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Start with SH PRS
load("./icd10_level_2_chapter_1-17_vs_sh_prs_test_nonrelatives_var1.RData")
load("./opcs_level_2_no_xyz_vs_sh_prs_test_nonrelatives_var1.RData")
ls()

thr <- 0.05/(3*(length(sh_test_nonr_var1_cov) + length(opcs_sh_test_nonr_var1) + 6)) # p = 5.50055e-05 is a significance threshold for glm() analyses; 6 - is a number of pain traits, glm results for them were included later

# ICD10
pval_icd <- lapply(sh_test_nonr_var1_cov, function(x) x[, 4])
pval_icd_i <- which(pval_icd < thr)

or_icd <- lapply(sh_test_nonr_var1_cov, function(x) exp(x[, 1]))
summary(unlist(or_icd))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.8915  1.0815  1.1253  1.1343  1.1769  1.3426
beta_icd <- lapply(sh_test_nonr_var1_cov, function(x) abs(x[, 1]))
or_icd_i <- which(beta_icd > log(2))

sh_icd_i <- intersect(pval_icd_i, or_icd_i) # 0 elemets
#sh_icd <- bp_test_nonr_var1_cov[bp_icd_i] # glm() results for ICD10 codes against SH PRS, significant and abs(beta) > log(2)
#or_sh_icd <- or_icd[sh_icd_i]
sh_icd <- sh_test_nonr_var1_cov[pval_icd_i] # glm() results for ICD10 codes against SH PRS, significant; N = 92
or_sh_icd <- or_icd[pval_icd_i]

readme_sh_icd <- "glm() results for ICD10 codes level 2, chapter I - XVII only, against SH PRS, significant"
save(readme_sh_icd, or_sh_icd, sh_icd, file = "./glm_icd_level_2_chapter_1-17_sh_prs_filtered_test_nonrelatives.RData")

# OPCS
pval_opcs <- lapply(opcs_sh_test_nonr_var1, function(x) x[, 4])
pval_opcs_i <- which(pval_opcs < thr)

or_opcs <- lapply(opcs_sh_test_nonr_var1, function(x) exp(x[, 1]))
summary(unlist(or_opcs))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9111  1.0471  1.1081  1.1102  1.1570  1.3683
beta_opcs <- lapply(opcs_sh_test_nonr_var1, function(x) abs(x[, 1]))
or_opcs_i <- which(beta_opcs > log(2))

sh_opcs_i <- intersect(pval_opcs_i, or_opcs_i) # 0 elements
#sh_opcs <- opcs_sh_test_nonr_var1[sh_opcs_i] # glm() results for ICD10 codes against SH PRS, significant and abs(beta) > log(2)
#or_sh_opcs <- or_opcs[sh_opcs_i]

sh_opcs <- opcs_sh_test_nonr_var1[pval_opcs_i] # glm() results for ICD10 codes against SH PRS, significant; N = 57
or_sh_opcs <- or_opcs[pval_opcs_i]

readme_sh_opcs <- "glm() results for OPCS codes level 2 (X, Y, Z OPCS codes excluded) against SH PRS, significant"
save(readme_sh_opcs, or_sh_opcs, sh_opcs, file = "./glm_opcs_level_2_no_xyz_sh_prs_filtered_test_nonrelatives.RData")



