# Aim of this script is to filter glm() results for ICD10 and OPCS codes
# against SH PRS binary coded

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Start with BP PRS
load("./icd10_vs_sh_prs.RData")
load("./opcs_vs_sh_prs.RData")
ls()

thr <- 0.05/(3*(length(sh_bin_cov) + length(opcs_sh_bin))) # p = 3.646973e-05 is a significance threshold for glm() analyses

# ICD10
pval_icd <- lapply(sh_bin_cov, function(x) x[, 4])
pval_icd_i <- which(pval_icd < thr)

or_icd <- lapply(sh_bin_cov, function(x) exp(x[, 1]))
beta_icd <- lapply(sh_bin_cov, function(x) abs(x[, 1]))
or_icd_i <- which(beta_icd > log(2))

sh_icd_i <- intersect(pval_icd_i, or_icd_i)
sh_icd <- sh_bin_cov[sh_icd_i] # glm() results for ICD10 codes against SH PRS, significant and abs(beta) > log(2)
or_sh_icd <- or_icd[sh_icd_i]

readme_sh_icd <- "glm() results for ICD10 codes against SH PRS, significant and abs(beta) > log(2)"
save(readme_sh_icd, or_sh_icd, sh_icd, file = "./glm_icd_sh_prs_filtered.RData")

# OPCS
pval_opcs <- lapply(opcs_sh_bin, function(x) x[, 4])
pval_opcs_i <- which(pval_opcs < thr)

or_opcs <- lapply(opcs_sh_bin, function(x) exp(x[, 1]))
beta_opcs <- lapply(opcs_sh_bin, function(x) abs(x[, 1]))
or_opcs_i <- which(beta_opcs > log(2))

sh_opcs_i <- intersect(pval_opcs_i, or_opcs_i)
sh_opcs <- opcs_sh_bin[sh_opcs_i] # glm() results for ICD10 codes against BP PRS, significant and abs(beta) > log(2)
or_sh_opcs <- or_opcs[sh_opcs_i]

readme_sh_opcs <- "glm() results for OPCS codes against SH PRS, significant and abs(beta) > log(2)"
save(readme_sh_opcs, or_sh_opcs, sh_opcs, file = "./glm_opcs_sh_prs_filtered.RData")



