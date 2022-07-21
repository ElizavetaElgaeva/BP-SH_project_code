# Aim of this script is to filter the glm() results
# for binary coded PRS traits against ICD10 and OPCS4 combined to the level 2
# Note: whole sample, cases, non-relatieves only; ICD10 codes from chapters I-XVII only; no X, Y, Z OPCS codes

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load ICD10 glm results
load("icd10_level_2_chapter1-17_vs_bin_prs_10_sh_all_cases.RData")
load("icd10_level_2_chapter1-17_vs_bin_prs_10_bp_sh_all_cases.RData")
load("icd10_level_2_chapter1-17_vs_bin_prs_90_bp_sh_all_cases.RData")
load("icd10_level_2_chapter1-17_vs_bin_sh_prs_90_sh_all_cases.RData")

# Load OPCS glm results
load("opcs_level_2_no_xyz_vs_bin_prs_10_sh_all_cases.RData")
load("opcs_level_2_no_xyz_bin_prs_10_bp_sh_all_cases.RData")
load("opcs_level_2_no_xyz_vs_bin_prs_90_sh_all_cases.RData")
load("opcs_level_2_no_xyz_vs_bin_prs_90_bp_sh_all_cases.RData")

ls()

# Set significance threshold
thr <- 0.05 / (4*(length(if_10_sh_vs_icd_l2) + length(if_90_sh_vs_opcs_l2))) # 3.434066e-05; length(if_10_sh_vs_icd_l2) = 199, length(if_90_sh_vs_opcs_l2) = 165 (all nonrel cases)

# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.1 quantile of SH PRS 
pval_icd_10_sh <- lapply(if_10_sh_vs_icd_l2, function(x) x[, 4])
pval_icd_10_sh_i <- which(pval_icd_10_sh < thr) # B96 E03 E11 E66 E78 F17 F32 G56 I10 I20 I21 I25 I73 J44 J45 K20 K21 K29 K44 K52 K57 K62 M13 M15 M17 M19 M23 M25 M47 M54 M75 M79
summary(unlist(pval_icd_10_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000000 0.0009134 0.0498205 0.1917959 0.3206108 0.9375529
icd_l2_vs_10_sh_sign <- do.call(rbind.data.frame, if_10_sh_vs_icd_l2[pval_icd_10_sh_i])
readme_icd_l2_vs_10_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, level 2, chapters 1-17, as dependent variable from binary coded SH PRS (10th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(icd_l2_vs_10_sh_sign, readme_icd_l2_vs_10_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_level_2_ch1-17_vs_10_sh_all_nonr_cases_signif.RData")

or_icd_10_sh <- lapply(if_10_sh_vs_icd_l2, function(x) exp(x[, 1]))
summary(unlist(or_icd_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4804  0.6979  0.7861  0.8079  0.8950  1.5171
beta_icd_10_sh <- lapply(if_10_sh_vs_icd_l2, function(x) abs(x[, 1]))
or_icd_10_sh_i <- which(beta_icd_10_sh > log(2)) 
intersect(pval_icd_10_sh_i, or_icd_10_sh_i) # 0


# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.1 quantile of BP-SH PRS
pval_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd_l2, function(x) x[, 4])
pval_icd_10_bp_sh_i <- which(pval_icd_10_bp_sh < thr) # E14 I10 I20 J45 K21 M13 M17 M19 M25
summary(unlist(pval_icd_10_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.02986 0.20235 0.33375 0.64862 0.99990
icd_l2_vs_10_bp_sh_sign <- do.call(rbind.data.frame, if_10_bp_sh_vs_icd_l2[pval_icd_10_bp_sh_i])
readme_icd_l2_vs_10_bp_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, level 2, chapters 1-17, as dependent variable from binary coded BP-SH PRS (10th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(icd_l2_vs_10_bp_sh_sign, readme_icd_l2_vs_10_bp_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_level_2_ch1-17_vs_10_bp_sh_all_nonr_cases_signif.RData")

or_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd_l2, function(x) exp(x[, 1]))
summary(unlist(or_icd_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7975  1.0126  1.1301  1.1213  1.2206  1.6414
beta_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd_l2, function(x) abs(x[, 1]))
or_icd_10_bp_sh_i <- which(beta_icd_10_bp_sh > log(2))
intersect(pval_icd_10_bp_sh_i, or_icd_10_bp_sh_i) # 0


# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.9 quantile of SH PRS
pval_icd_90_sh <- lapply(if_90_sh_vs_icd_l2, function(x) x[, 4])
pval_icd_90_sh_i <- which(pval_icd_90_sh < thr) # A09 B96 D17 E03 E11 E66 E78 F17 F32 F41 G40 G43 G47 G56 I10 I20 I25 I73 I84 I95 J22 J34 J43 J44 J45 K20 K21 K26 K29 K30 K31 K42 K44 K51 K52 K57 K58 K59 K62 K63 K76 K80 K85 K92 L03 M06 M13 M15 M16 M17 M19 M23 M25 M35 M47 M48 M50 M51 M54 M70 M75 M77 M79 N18 N28 N39 N92
summary(unlist(pval_icd_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000000 0.0000012 0.0016496 0.1073396 0.0802277 0.9351126
icd_l2_vs_90_sh_sign <- do.call(rbind.data.frame, if_90_sh_vs_icd_l2[pval_icd_90_sh_i])
readme_icd_l2_vs_90_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, level 2, chapters 1-17, as dependent variable from binary coded SH PRS (90th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(icd_l2_vs_90_sh_sign, readme_icd_l2_vs_90_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_level_2_ch1-17_vs_90_sh_all_nonr_cases_signif.RData")

or_icd_90_sh <- lapply(if_90_sh_vs_icd_l2, function(x) exp(x[, 1]))
summary(unlist(or_icd_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.6646  1.1878  1.4009  1.3713  1.5294  1.9088
beta_icd_90_sh <- lapply(if_90_sh_vs_icd_l2, function(x) abs(x[, 1]))
or_icd_90_sh_i <- which(beta_icd_90_sh > log(2))
intersect(pval_icd_90_sh_i, or_icd_90_sh_i) # 0


# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.9 quantile of BP-SH PRS
pval_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd_l2, function(x) x[, 4])
pval_icd_90_bp_sh_i <- which(pval_icd_90_bp_sh < thr) # E66 K21 M13 M16 M17 M19 M23 M25 M47 M75
summary(unlist(pval_icd_90_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.03094 0.17861 0.30417 0.53754 0.99837
icd_l2_vs_90_bp_sh_sign <- do.call(rbind.data.frame, if_90_bp_sh_vs_icd_l2[pval_icd_90_bp_sh_i])
readme_icd_l2_vs_90_bp_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, level 2, chapters 1-17, as dependent variable from binary coded BP-SH PRS (90th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(icd_l2_vs_90_bp_sh_sign, readme_icd_l2_vs_90_bp_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_level_2_ch1-17_vs_90_bp_sh_all_nonr_cases_signif.RData")

or_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd_l2, function(x) exp(x[, 1]))
summary(unlist(or_icd_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4516  0.7820  0.8708  0.8759  0.9633  1.4118
beta_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd_l2, function(x) abs(x[, 1]))
or_icd_90_bp_sh_i <- which(beta_icd_90_bp_sh > log(2))
intersect(pval_icd_90_bp_sh_i, or_icd_90_bp_sh_i) # 0


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.1 quantile of SH PRS
pval_opcs_10_sh <- lapply(if_10_sh_vs_opcs_l2, function(x) x[, 4])
pval_opcs_10_sh_i <- which(pval_opcs_10_sh < thr) # A65 G45 H22 H25 K63 M45 V54 W40 W82 W83 W84 W90
summary(unlist(pval_opcs_10_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.00756 0.10029 0.23164 0.38258 0.96367
opcs_l2_vs_10_sh_sign <- do.call(rbind.data.frame, if_10_sh_vs_opcs_l2[pval_opcs_10_sh_i])
readme_opcs_l2_vs_10_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, level 2, chapters X, Y, Z excluded, as dependent variable from binary coded SH PRS (10th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(opcs_l2_vs_10_sh_sign, readme_opcs_l2_vs_10_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_level_2_no_xyz_vs_10_sh_all_nonr_cases_signif.RData")

or_opcs_10_sh <- lapply(if_10_sh_vs_opcs_l2, function(x) exp(x[, 1]))
summary(unlist(or_opcs_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4854  0.7066  0.8124  0.8268  0.9302  1.4019
beta_opcs_10_sh <- lapply(if_10_sh_vs_opcs_l2, function(x) abs(x[, 1]))
or_opcs_10_sh_i <- which(beta_opcs_10_sh > log(2))
intersect(or_opcs_10_sh_i, pval_opcs_10_sh_i) # W83


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.1 quantile of BP-SH PRS
pval_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs_l2, function(x) x[, 4])
pval_opcs_10_bp_sh_i <- which(pval_opcs_10_bp_sh < thr) # E25 G45 W40 W82
summary(unlist(pval_opcs_10_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.08335 0.39135 0.40405 0.66363 0.99438
opcs_l2_vs_10_bp_sh_sign <- do.call(rbind.data.frame, if_10_bp_sh_vs_opcs_l2[pval_opcs_10_bp_sh_i])
readme_opcs_l2_vs_10_bp_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, level 2, chapters X, Y, Z excluded, as dependent variable from binary coded BP-SH PRS (10th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(opcs_l2_vs_10_bp_sh_sign, readme_opcs_l2_vs_10_bp_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_level_2_no_xyz_vs_10_bp_sh_all_nonr_cases_signif.RData")

or_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs_l2, function(x) exp(x[, 1]))
summary(unlist(or_opcs_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.6562  0.9845  1.0630  1.0908  1.1893  1.5472
beta_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs_l2, function(x) abs(x[, 1]))
or_opcs_10_bp_sh_i <- which(beta_opcs_10_bp_sh > log(2))
intersect(pval_opcs_10_bp_sh_i,or_opcs_10_bp_sh_i) # 0


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.9 quantile of SH PRS
pval_opcs_90_sh <- lapply(if_90_sh_vs_opcs_l2, function(x) x[, 4])
pval_opcs_90_sh_i <- which(pval_opcs_90_sh < thr) # A52 A57 A65 A73 E25 E36 G19 G45 H22 H25 J18 K63 M45 O16 O29 Q07 T24 T62 T79 U05 U21 V29 V54 V55 W06 W37 W40 W80 W82 W83 W84 W85 W87 W90 W91
summary(unlist(pval_opcs_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000000 0.0000918 0.0449982 0.2308736 0.4382607 0.9844360
opcs_l2_vs_90_sh_sign <- do.call(rbind.data.frame, if_90_sh_vs_opcs_l2[pval_opcs_90_sh_i])
readme_opcs_l2_vs_90_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, level 2, chapters X, Y, Z excluded, as dependent variable from binary coded SH PRS (90th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(opcs_l2_vs_90_sh_sign, readme_opcs_l2_vs_90_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_level_2_no_xyz_vs_90_sh_all_nonr_cases_signif.RData")

or_opcs_90_sh <- lapply(if_90_sh_vs_opcs_l2, function(x) exp(x[, 1]))
summary(unlist(or_opcs_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7013  1.0895  1.2693  1.2830  1.4158  2.0842
beta_opcs_90_sh <- lapply(if_90_sh_vs_opcs_l2, function(x) abs(x[, 1]))
or_opcs_90_sh_i <- which(beta_opcs_90_sh > log(2))
intersect(pval_opcs_90_sh_i, or_opcs_90_sh_i) # E25 T62 V29


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.9 quantile both of BP-SH PRS
pval_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs_l2, function(x) x[, 4])
pval_opcs_90_bp_sh_i <- which(pval_opcs_90_bp_sh < thr) # G45 H22 W82 W40
summary(unlist(pval_opcs_90_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.06423 0.33666 0.37878 0.64661 0.99492
opcs_l2_vs_90_bp_sh_sign <- do.call(rbind.data.frame, if_90_bp_sh_vs_opcs_l2[pval_opcs_90_bp_sh_i])
readme_opcs_l2_vs_90_bp_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, level 2, chapters X, Y, Z excluded, as dependent variable from binary coded BP-SH PRS (90th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(opcs_l2_vs_90_bp_sh_sign, readme_opcs_l2_vs_90_bp_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_level_2_no_xyz_vs_90_bp_sh_all_nonr_cases_signif.RData")

or_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs_l2, function(x) exp(x[, 1]))
summary(unlist(or_opcs_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.5202  0.8079  0.9206  0.9169  0.9991  1.5867
beta_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs_l2, function(x) abs(x[, 1]))
or_opcs_90_bp_sh_i <- which(beta_opcs_90_bp_sh > log(2))
intersect(pval_opcs_90_bp_sh_i, or_opcs_90_bp_sh_i) # 0


