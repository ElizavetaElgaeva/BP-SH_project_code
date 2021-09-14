# Aim of this script is to filter glm() results for ICD10 and OPCS codes
# against BP PRS binary coded

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Start with BP PRS
load("./icd10_vs_bp_prs.RData")
load("./opcs_vs_bp_prs.RData")
ls()

thr <- 0.05/(3*(length(bp_bin_cov) + length(opcs_bp_bin))) # p = 3.646973e-05 is a significance threshold for glm() analyses

# ICD10
pval_icd <- lapply(bp_bin_cov, function(x) x[, 4])
pval_icd_i <- which(pval_icd < thr)

or_icd <- lapply(bp_bin_cov, function(x) exp(x[, 1]))
beta_icd <- lapply(bp_bin_cov, function(x) abs(x[, 1]))
or_icd_i <- which(beta_icd > log(2))

bp_icd_i <- intersect(pval_icd_i, or_icd_i)
bp_icd <- bp_bin_cov[bp_icd_i] # glm() results for ICD10 codes against BP PRS, significant and abs(beta) > log(2)
or_bp_icd <- or_icd[bp_icd_i]

readme_bp_icd <- "glm() results for ICD10 codes against BP PRS, significant and abs(beta) > log(2)"
save(readme_bp_icd, or_bp_icd, bp_icd, file = "./glm_icd_bp_prs_filtered.RData")

# OPCS
pval_opcs <- lapply(opcs_bp_bin, function(x) x[, 4])
pval_opcs_i <- which(pval_opcs < thr)

or_opcs <- lapply(opcs_bp_bin, function(x) exp(x[, 1]))
beta_opcs <- lapply(opcs_bp_bin, function(x) abs(x[, 1]))
or_opcs_i <- which(beta_opcs > log(2))

bp_opcs_i <- intersect(pval_opcs_i, or_opcs_i)
bp_opcs <- opcs_bp_bin[bp_opcs_i] # glm() results for ICD10 codes against BP PRS, significant and abs(beta) > log(2)
or_bp_opcs <- or_opcs[bp_opcs_i]

readme_bp_opcs <- "glm() results for OPCS codes against BP PRS, significant and abs(beta) > log(2)"
save(readme_bp_opcs, or_bp_opcs, bp_opcs, file = "./glm_opcs_bp_prs_filtered.RData")



