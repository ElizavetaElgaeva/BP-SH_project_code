# Aim of this script is to filter glm() results for ICD10 and OPCS codes
# against BP-SH PRS binary coded

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Start with BP PRS
load("./icd10_vs_bp_sh_prs.RData")
load("./opcs_vs_bp_sh_prs.RData")
ls()

thr <- 0.05/(3*(length(bp_sh_bin_cov) + length(opcs_bp_sh_bin))) # p = 3.646973e-05 is a significance threshold for glm() analyses

# ICD10
pval_icd <- lapply(bp_sh_bin_cov, function(x) x[, 4])
pval_icd_i <- which(pval_icd < thr)

or_icd <- lapply(bp_sh_bin_cov, function(x) exp(x[, 1]))
beta_icd <- lapply(bp_sh_bin_cov, function(x) abs(x[, 1]))
or_icd_i <- which(beta_icd > log(2))

bp_sh_icd_i <- intersect(pval_icd_i, or_icd_i) # empty

#bp_sh_icd <- bp_sh_bin_cov[bp_sh_icd_i] # glm() results for ICD10 codes against BP PRS, significant and OR > 2
#or_bp_sh_icd <- or_icd[bp_sh_icd_i]

#readme_bp_sh_icd <- "glm() results for ICD10 codes against BP-SH PRS, significant and OR > 2"
#save(readme_bp_sh_icd, or_bp_sh_icd, bp_sh_icd, file = "./glm_icd_bp_sh_prs_filtered.RData")

# OPCS
pval_opcs <- lapply(opcs_bp_sh_bin, function(x) x[, 4])
pval_opcs_i <- which(pval_opcs < thr)

or_opcs <- lapply(opcs_bp_sh_bin, function(x) exp(x[, 1]))
beta_opcs <- lapply(opcs_bp_sh_bin, function(x) abs(x[, 1]))
or_opcs_i <- which(beta_opcs > log(2))

bp_sh_opcs_i <- intersect(pval_opcs_i, or_opcs_i) # empty


## Select top-5 abs(beta) values among significant glm() results
# ICD10
bp_sh_icd_sig <- bp_sh_bin_cov[pval_icd_i]
beta_icd_sig <- lapply(bp_sh_icd_sig, function(x) abs(x[, 1]))
summary(unlist(beta_icd_sig))
icd_beta_top_sig <- names(sort(unlist(beta_icd_sig), decreasing = T)[1:5])
bp_sh_icd_beta <- sort(unlist(beta_icd_sig), decreasing = T)[1:5]
bp_sh_icd_beta_top_sig <- bp_sh_icd_sig[icd_beta_top_sig]
bp_sh_icd_or <- lapply(bp_sh_icd_beta_top_sig, function(x) exp(x[, 1]))
readme_bp_sh_icd <- "glm() results for ICD10 codes against BP-SH PRS, top-5 abs(beta) values among significant results"
save(readme_bp_sh_icd, bp_sh_icd_beta_top_sig, bp_sh_icd_or, file = "./glm_icd_bp_sh_prs_filtered.RData")


bp_sh_opcs_sig <- opcs_bp_sh_bin[pval_opcs_i]
beta_opcs_sig <- lapply(bp_sh_opcs_sig, function(x) abs(x[, 1]))
summary(unlist(beta_opcs_sig))
opcs_beta_top_sig <- names(sort(unlist(beta_opcs_sig), decreasing = T)[1:5])
bp_sh_opcs_beta <- sort(unlist(beta_opcs_sig), decreasing = T)[1:5]
bp_sh_opcs_beta_top_sig <- bp_sh_opcs_sig[opcs_beta_top_sig]
bp_sh_opcs_or <- lapply(bp_sh_opcs_beta_top_sig, function(x) exp(x[, 1]))
readme_bp_sh_opcs <- "glm() results for OPCS codes against BP-SH PRS, top-5 abs(beta) values among significant results"
save(readme_bp_sh_opcs, bp_sh_opcs_beta_top_sig, bp_sh_opcs_or, file = "./glm_opcs_bp_sh_prs_filtered.RData")

# ICD10 & OPCS joint results
c(bp_sh_icd_or, bp_sh_opcs_or)
sort(c(unlist(bp_sh_icd_or), unlist(bp_sh_opcs_or)), decreasing = T)[1:5] # top-5 among ICD10 and OPCS
icd_opcs_beta_top_sig <- names(sort(c(unlist(bp_sh_icd_or), unlist(bp_sh_opcs_or)), decreasing = T)[1:5])
icd_opcs <- c(bp_sh_opcs_beta_top_sig, bp_sh_icd_beta_top_sig)
icd_opcs <- icd_opcs[icd_opcs_beta_top_sig]
icd_opcs_or <- lapply(icd_opcs, function(x) exp(x[, 1]))
readme_icd_opcs <- "glm() results for ICD10 and OPCS codes against BP-SH PRS, top-5 abs(beta) values among significant results"
save(readme_icd_opcs, icd_opcs, icd_opcs_or, file = "./glm_icd_opcs_bp_sh_prs_filtered.RData")


#bp_sh_opcs <- opcs_bp_sh_bin[bp_sh_opcs_i] # glm() results for ICD10 codes against BP-SH PRS, significant and OR > 2
#or_bp_sh_opcs <- or_opcs[bp_sh_opcs_i]

#readme_bp_sh_opcs <- "glm() results for OPCS codes against BP-SH PRS, significant and OR > 2"
#save(readme_bp_sh_opcs, or_bp_sh_opcs, bp_sh_opcs, file = "./glm_opcs_bp_sh_prs_filtered.RData")



