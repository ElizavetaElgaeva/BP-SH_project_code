# Aim of this script is to filter the glm() results 
# for binary coded PRS traits (non-overlapping 0.1 and0.9 quantiles of SH and BP-SH PRS) against ICD10 and OPCS4 not-combined
# Note: test sample, cases, non-relatieves only; ICD10 codes from chapters I-XVII only; no X, Y, Z OPCS codes

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load ICD10 glm results
load("icd10_chapter1-17vs_bin_sh_prs_90_sh.RData")
load("icd10_chapter1-17_vs_bin_prs_10_bp_sh.RData")
load("icd10_chapter1-17_vs_bin_prs_10_sh.RData")
load("icd10_chapter1-17_vs_bin_prs_90_bp_sh.RData")
#load("icd10_vs_bin_sh_prs_90_sh.RData") # non-overlappping quantiles, test sample, not-combined but all codes
#load("icd10_level_2_chapter_1-17_vs_bin_sh_prs_90_sh.RData") # non-overlappping quantiles, test sample, combined codes till chapter R
 
# Load OPCS glm results
load("opcs_no_xyz_vs_bin_prs_90_bp_sh.RData")
load("opcs_no_xyz_bin_prs_10_bp_sh.RData")
load("opcs_no_xyz_vs_bin_prs_90_sh.RData")
load("opcs_no_xyz_vs_bin_prs_10_sh.RData")
#load("opcs_vs_bin_prs_10_sh.RData") # non-overlappping quantiles, test sample, not-combined but all codes
#load("opcs_bin_prs_10_bp_sh.RData")

ls()

# Set significance threshold
thr <- 0.05 / (4*(length(if_10_sh_vs_icd) + length(if_90_sh_vs_opcs))) # 4.464286e-05; length(if_10_sh_vs_icd) = 170, length(if_90_bp_sh_vs_opcs) = 110

# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.1 quantile of SH PRS 
pval_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) x[, 4])
pval_icd_10_sh_i <- which(pval_icd_10_sh < thr) # I10
summary(unlist(pval_icd_10_sh))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0000013 0.0579529 0.2621130 0.3374100 0.5573738 0.9877788
icd_vs_10_sh_sign <- do.call(rbind.data.frame, if_10_sh_vs_icd[pval_icd_10_sh_i])
readme_icd_vs_10_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, not-combined, chapters 1-17, as dependent variable from binary coded SH PRS (10th percentile - 1, other - 0), test sample of nonrelatives, cases"
save(icd_vs_10_sh_sign, readme_icd_vs_10_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_ch1-17_vs_10_sh_test_nonr_cases_signif.RData")

or_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3203  0.7025  0.8148  0.8339  0.9326  1.5842
beta_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) abs(x[, 1]))
or_icd_10_sh_i <- which(beta_icd_10_sh > log(2)) # J348  K29 I517
intersect(pval_icd_10_sh_i, or_icd_10_sh_i) # 0 elements


# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.1 quantile of BP-SH PRS
pval_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) x[, 4])
pval_icd_10_bp_sh_i <- which(pval_icd_10_bp_sh < thr) # 0 elements
summary(unlist(pval_icd_10_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0005579 0.2160958 0.5131061 0.4774730 0.7023528 0.9844495

or_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3525  0.8676  1.0047  1.0017  1.1099  1.8075
beta_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) abs(x[, 1]))
or_icd_10_bp_sh_i <- which(beta_icd_10_bp_sh > log(2))


# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.9 quantile of SH PRS
pval_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) x[, 4])
pval_icd_90_sh_i <- which(pval_icd_90_sh < thr) # K219 M754 M545 M179 K449 F171
summary(unlist(pval_icd_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000004 0.0439470 0.3050774 0.3417980 0.5706464 0.9807432
icd_vs_90_sh_sign <- do.call(rbind.data.frame, if_90_sh_vs_icd[pval_icd_90_sh_i])
readme_icd_vs_90_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, not-combined, chapters 1-17, as dependent variable from binary coded SH PRS (90th percentile - 1, other - 0), test sample of nonrelatives, cases"
save(icd_vs_90_sh_sign, readme_icd_vs_90_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_ch1-17_vs_90_sh_test_nonr_cases_signif.RData")

or_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.652   1.009   1.167   1.181   1.351   1.919
beta_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) abs(x[, 1]))
or_icd_90_sh_i <- which(beta_icd_90_sh > log(2))
intersect(pval_icd_90_sh_i, or_icd_90_sh_i) # 0


# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.9 quantile of BP-SH PRS
pval_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) x[, 4])
pval_icd_90_bp_sh_i <- which(pval_icd_90_bp_sh < thr) # 0
summary(unlist(pval_icd_90_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.004072 0.284498 0.508252 0.521366 0.776025 0.990672

or_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4502  0.8938  1.0031  1.0015  1.1128  1.7410
beta_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) abs(x[, 1]))
or_icd_90_bp_sh_i <- which(beta_icd_90_bp_sh > log(2))


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.1 quantile of SH PRS
pval_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) x[, 4])
pval_opcs_10_sh_i <- which(pval_opcs_10_sh < thr) # G451
summary(unlist(pval_opcs_10_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000164 0.1272975 0.3662814 0.3917125 0.5751957 0.9749227
opcs_vs_10_sh_sign <- do.call(rbind.data.frame, if_10_sh_vs_opcs[pval_opcs_10_sh_i])
readme_opcs_vs_10_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, not-combined, chapters X, Y, Z excluded, as dependent variable from binary coded SH PRS (10th percentile - 1, other - 0), test sample of nonrelatives, cases"
save(opcs_vs_10_sh_sign, readme_opcs_vs_10_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_10_sh_test_nonr_cases_signif.RData")

or_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) exp(x[, 1]))
summary(unlist(or_opcs_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4816  0.7676  0.8585  0.8871  0.9469  2.9243
beta_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) abs(x[, 1]))
or_opcs_10_sh_i <- which(beta_opcs_10_sh > log(2))
intersect(or_opcs_10_sh_i, pval_opcs_10_sh_i) # 0


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.1 quantile of BP-SH PRS
pval_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) x[, 4])
pval_opcs_10_bp_sh_i <- which(pval_opcs_10_bp_sh < thr) # 0 elements
summary(unlist(pval_opcs_10_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00898 0.17382 0.47021 0.46487 0.75514 0.97118

or_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) exp(x[, 1]))
summary(unlist(or_opcs_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3244  0.8230  0.9464  0.9325  1.0487  1.5693
beta_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) abs(x[, 1]))
or_opcs_10_bp_sh_i <- which(beta_opcs_10_bp_sh > log(2))


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.9 quantile of SH PRS
pval_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) x[, 4])
pval_opcs_90_sh_i <- which(pval_opcs_90_sh < thr) # G459 G451 H229 W401 
summary(unlist(pval_opcs_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000  0.1062  0.3175  0.3988  0.6731  0.9880
opcs_vs_90_sh_sign <- do.call(rbind.data.frame, if_90_sh_vs_opcs[pval_opcs_90_sh_i])
readme_opcs_vs_90_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, not-combined, chapters X, Y, Z excluded, as dependent variable from binary coded SH PRS (90th percentile - 1, other - 0), test sample of nonrelatives, cases"
save(opcs_vs_90_sh_sign, readme_opcs_vs_90_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_90_sh_test_nonr_cases_signif.RData")

or_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) exp(x[, 1]))
summary(unlist(or_opcs_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.3534  0.9920  1.1162  1.1436  1.3368  1.9773
beta_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) abs(x[, 1]))
or_opcs_90_sh_i <- which(beta_opcs_90_sh > log(2))
intersect(pval_opcs_90_sh_i, or_opcs_90_sh_i) # 0


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.9 quantile of BP-SH PRS
pval_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) x[, 4])
pval_opcs_90_bp_sh_i <- which(pval_opcs_90_bp_sh < thr) # 0 elements
summary(unlist(pval_opcs_90_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.004296 0.265883 0.537836 0.527552 0.798690 0.998465

or_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) exp(x[, 1]))
summary(unlist(or_opcs_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2729  0.8943  0.9929  0.9979  1.1041  1.6843
beta_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) abs(x[, 1]))
or_opcs_90_bp_sh_i <- which(beta_opcs_90_bp_sh > log(2))

