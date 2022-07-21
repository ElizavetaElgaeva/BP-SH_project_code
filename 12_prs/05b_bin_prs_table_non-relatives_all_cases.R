# Aim of this script is to create a table with glm() results for binary coded SH/BP-SH PRS (non-overlapping 0.1 and 0.9 quantiles)
# (Whole sample, cases, non-relatieves only; standardized PRS; ICD10 and OPCS codes combined to the level 2)
# Note: ICD10 from chapter I - XVII only; X, Y, Z OPCS codes excluded

library(data.table)
library(dplyr)

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load full glm() results
# SH
load("icd10_level_2_chapter1-17_vs_bin_prs_10_sh_all_cases.RData")
load("icd10_level_2_chapter1-17_vs_bin_sh_prs_90_sh_all_cases.RData")

load("opcs_level_2_no_xyz_vs_bin_prs_10_sh_all_cases.RData")
load("opcs_level_2_no_xyz_vs_bin_prs_90_sh_all_cases.RData")

# BP-SH
load("icd10_level_2_chapter1-17_vs_bin_prs_10_bp_sh_all_cases.RData")
load("icd10_level_2_chapter1-17_vs_bin_prs_90_bp_sh_all_cases.RData")

load("opcs_level_2_no_xyz_bin_prs_10_bp_sh_all_cases.RData")
load("opcs_level_2_no_xyz_vs_bin_prs_90_bp_sh_all_cases.RData")


# Load glm results filtered by statistical significance
# SH PRS
load("icd_level_2_ch1-17_vs_10_sh_all_nonr_cases_signif.RData")
load("icd_level_2_ch1-17_vs_90_sh_all_nonr_cases_signif.RData")

load("opcs_level_2_no_xyz_vs_10_sh_all_nonr_cases_signif.RData")
load("opcs_level_2_no_xyz_vs_90_sh_all_nonr_cases_signif.RData")

# BP-SH PRS
load("icd_level_2_ch1-17_vs_10_bp_sh_all_nonr_cases_signif.RData")
load("icd_level_2_ch1-17_vs_90_bp_sh_all_nonr_cases_signif.RData")

load("opcs_level_2_no_xyz_vs_10_bp_sh_all_nonr_cases_signif.RData")
load("opcs_level_2_no_xyz_vs_90_bp_sh_all_nonr_cases_signif.RData")

ls()

# Select ICD10 and OPCS codes of interest
codes <- rownames(icd_l2_vs_10_sh_sign)
codes <- c(codes, rownames(icd_l2_vs_90_sh_sign))
codes <- c(codes, rownames(icd_l2_vs_10_bp_sh_sign))
codes <- c(codes, rownames(icd_l2_vs_90_bp_sh_sign))
codes <- unique(codes)

codes2 <- rownames(opcs_l2_vs_10_sh_sign)
codes2 <- c(codes2, rownames(opcs_l2_vs_90_sh_sign))
codes2 <- c(codes2, rownames(opcs_l2_vs_10_bp_sh_sign))
codes2 <- c(codes2, rownames(opcs_l2_vs_90_bp_sh_sign))
codes2 <- unique(codes2)

intersect(codes, codes2) # K63
all_codes <- c(codes, codes2) # length(codes) = 69, length(codes2) = 35
icd_ind <- c(1:length(codes))

# Create a table
tab <- as.data.frame(matrix(ncol = 20, nrow = length(all_codes)))
colnames(tab) <- c("nomenclature", "code", "description", "10_sh_or", "10_sh_beta", "10_sh_se", "10_sh_p", "90_sh_or", "90_sh_beta", "90_sh_se", "90_sh_p", "10_bp_sh_or", "10_bp_sh_beta", "10_bp_sh_se", "10_bp_sh_p", "90_bp_sh_or", "90_bp_sh_beta", "90_bp_sh_se", "90_bp_sh_p", "code_prevalence")
tab$code <- all_codes
tab$nomenclature[icd_ind] <- "ICD10"
tab$nomenclature[-icd_ind] <- "OPCS4"

# Add description
icd_des <- readxl::read_excel("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/ICD10_coding19.xlsx", sheet = 1)
opcs_des <- readxl::read_excel("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/OPCS_coding240.xlsx", sheet = 1)

i_icd_des <- match(codes, icd_des$"coding")
i_opcs_des <- match(codes2, opcs_des$"coding")

icd_des <- icd_des[i_icd_des, c("coding", "meaning")]
opcs_des <- opcs_des[i_opcs_des, c("coding", "meaning")]
icd_opcs_des <- rbind(icd_des, opcs_des)
table(tab$code == icd_opcs_des$"coding")
tab$description <- icd_opcs_des$"meaning"


sh_10_full <- c(if_10_sh_vs_icd_l2[codes], if_10_sh_vs_opcs_l2[codes2])
tab$"10_sh_or" <- sapply(sh_10_full, function(x) exp(x[, 1]))
tab$"10_sh_beta" <- sapply(sh_10_full, function(x) x[, 1])
tab$"10_sh_se" <- sapply(sh_10_full, function(x) x[, 2])
tab$"10_sh_p" <- sapply(sh_10_full, function(x) x[, 4])

sh_90_full <- c(if_90_sh_vs_icd_l2[codes], if_90_sh_vs_opcs_l2[codes2])
tab$"90_sh_or" <- sapply(sh_90_full, function(x) exp(x[, 1]))
tab$"90_sh_beta" <- sapply(sh_90_full, function(x) x[, 1])
tab$"90_sh_se" <- sapply(sh_90_full, function(x) x[, 2])
tab$"90_sh_p" <- sapply(sh_90_full, function(x) x[, 4])

bp_sh_10_full <- c(if_10_bp_sh_vs_icd_l2[codes], if_10_bp_sh_vs_opcs_l2[codes2])
tab$"10_bp_sh_or" <- sapply(bp_sh_10_full, function(x) exp(x[, 1]))
tab$"10_bp_sh_beta" <- sapply(bp_sh_10_full, function(x) x[, 1])
tab$"10_bp_sh_se" <- sapply(bp_sh_10_full, function(x) x[, 2])
tab$"10_bp_sh_p" <- sapply(bp_sh_10_full, function(x) x[, 4])

bp_sh_90_full <- c(if_90_bp_sh_vs_icd_l2[codes], if_90_bp_sh_vs_opcs_l2[codes2])
tab$"90_bp_sh_or" <- sapply(bp_sh_90_full, function(x) exp(x[, 1]))
tab$"90_bp_sh_beta" <- sapply(bp_sh_90_full, function(x) x[, 1])
tab$"90_bp_sh_se" <- sapply(bp_sh_90_full, function(x) x[, 2])
tab$"90_bp_sh_p" <- sapply(bp_sh_90_full, function(x) x[, 4])


load("//home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/prev_of_icd10_level_2_cases_all_nonrel_filtered.RData")
load("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/prev_opcs_level_2_cases_all_nonrel_filtered.RData")
prev_full <- c(prev_icd_l2_cases_f[codes], prev_opcs_l2_cases_f[codes2])
tab$"code_prevalence" <- prev_full

fwrite(tab, file = "glm_level_2_codes_bin_prs_results_preselected_codes_nonrel_all_cases.txt", sep = "\t", dec = ".")


