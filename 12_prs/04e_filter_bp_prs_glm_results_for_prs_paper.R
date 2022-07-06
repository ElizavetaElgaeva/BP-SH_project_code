# Aim of this script is to filter glm() results for ICD10 and OPCS codes (both combined to level 2)
# against standardized BP PRS (test sample non-relatieves only)
# Note: ICD10 from chapter I - XVII only; X, Y, Z OPCS codes excluded

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Start with BP PRS
load("./icd10_level_2_chapter_1-17_vs_bp_prs_test_nonrelatives_var1.RData")
load("./opcs_level_2_no_xyz_vs_bp_prs_test_nonrelatives_var1.RData")
ls()

thr <- 0.05/(length(bp_test_nonr_var1_cov) + length(opcs_bp_test_nonr_var1) + 6) # 0.0001650165

# ICD10
pval_icd <- lapply(bp_test_nonr_var1_cov, function(x) x[, 4])
pval_icd_i <- which(pval_icd < thr)

or_icd <- lapply(bp_test_nonr_var1_cov, function(x) exp(x[, 1]))
summary(unlist(or_icd))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.918   1.052   1.094   1.096   1.129   1.300
beta_icd <- lapply(bp_test_nonr_var1_cov, function(x) abs(x[, 1]))
or_icd_i <- which(beta_icd > log(2))

bp_icd_i <- intersect(pval_icd_i, or_icd_i) # 0 elemets
#bp_icd <- bp_test_nonr_var1_cov[bp_icd_i] # glm() results for ICD10 codes against BP PRS, significant and abs(beta) > log(2)
#or_bp_icd <- or_icd[bp_icd_i]
bp_icd <- bp_test_nonr_var1_cov[pval_icd_i] # glm() results for ICD10 codes against BP PRS, significant; N = 62
or_bp_icd <- or_icd[pval_icd_i]

readme_bp_icd <- "glm() results for ICD10 codes level 2, chapter I - XVII only, against BP PRS, significant, threshold for PRS paper"
save(readme_bp_icd, or_bp_icd, bp_icd, file = "./glm_icd_level_2_chapter_1-17_bp_prs_filtered_test_nonrelatives_for_PRS_paper.RData")

# OPCS
pval_opcs <- lapply(opcs_bp_test_nonr_var1, function(x) x[, 4])
pval_opcs_i <- which(pval_opcs < thr)

or_opcs <- lapply(opcs_bp_test_nonr_var1, function(x) exp(x[, 1]))
summary(unlist(or_opcs))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9269  1.0278  1.0776  1.0808  1.1183  1.3519
beta_opcs <- lapply(opcs_bp_test_nonr_var1, function(x) abs(x[, 1]))
or_opcs_i <- which(beta_opcs > log(2))

bp_opcs_i <- intersect(pval_opcs_i, or_opcs_i) # 0 elements
#bp_opcs <- opcs_bp_test_nonr_var1[bp_opcs_i] # glm() results for ICD10 codes against BP PRS, significant and abs(beta) > log(2)
#or_bp_opcs <- or_opcs[bp_opcs_i]

bp_opcs <- opcs_bp_test_nonr_var1[pval_opcs_i] # glm() results for ICD10 codes against BP PRS, significant; N = 37
or_bp_opcs <- or_opcs[pval_opcs_i]

readme_bp_opcs <- "glm() results for OPCS codes level 2 (X, Y, Z OPCS codes excluded) against BP PRS, significant, threshold for PRS paper"
save(readme_bp_opcs, or_bp_opcs, bp_opcs, file = "./glm_opcs_level_2_no_xyz_bp_prs_filtered_test_nonrelatives_for_PRS_paper.RData")



