# Aim of this script is to filter the glm() results 
# for binary coded PRS traits against ICD10 and OPCS4 level 2 codes
# Note: test sample, non-relatieves only; ICD10 codes from chapters I-XVII only; no X, Y, Z OPCS codes

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load ICD10 glm results
load("icd10_level_2_chapter_1-17_vs_bin_prs_10_10.RData")
load("icd10_level_2_chapter_1-17_vs_bin_prs_10_90.RData")
load("icd10_level_2_chapter_1-17_vs_bin_sh_prs_90_10.RData")
load("icd10_level_2_chapter_1-17_vs_bin_prs_90_90.RData")

#load("icd10_chapter1-17_vs_bin_prs_10_sh_all_cases.RData")

#load("icd10_chapter1-17vs_bin_sh_prs_90_sh.RData")
# load("icd10_vs_bin_sh_prs_90_sh.RData")
#load("/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter_1-17_vs_bin_sh_prs_90_sh.RData")
 
# Load OPCS glm results
load("opcs_level_2_no_xyz_vs_bin_prs_10_10.RData")
load("opcs_level_2_no_xyz_vs_bin_prs_10_90.RData")
load("opcs_level_2_no_xyz_vs_bin_sh_prs_90_10.RData")
load("opcs_level_2_no_xyz_vs_bin_prs_90_90.RData")

#load("opcs_no_xyz_vs_bin_prs_10_sh_all_cases.RData")

#load("opcs_no_xyz_vs_bin_prs_90_bp_sh.RData")
# load("opcs_vs_bin_prs_10_sh.RData")
#load("opcs_bin_prs_10_bp_sh.RData")
ls()

# Set significance threshold
thr <- 0.05 / (4*(length(if_10_10_vs_icd) + length(if_10_10_vs_opcs))) # 4.208754e-05; length(if_10_10_vs_icd) = 165, length(if_10_10_vs_opcs) = 132
#thr <- 0.05 / (4*(length(if_10_sh_vs_icd) + length(if_90_sh_vs_opcs))) # 4.464286e-05; length(if_10_sh_vs_icd) = 170, length(if_90_bp_sh_vs_opcs) = 110
#thr <- 0.05 / (4*(length(if_10_sh_vs_icd) + length(if_90_sh_vs_opcs))) # 3.140704e-05; length(if_10_sh_vs_icd) = 243, length(if_90_sh_vs_opcs) = 155 (all nonrel cases)

# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.1 quantile for both SH and BP-SH PRS 
pval_icd_10_10 <- lapply(if_10_10_vs_icd, function(x) x[, 4])
pval_icd_10_10_i <- which(pval_icd_10_10 < thr) # 0 elements
summary(unlist(pval_icd_10_10))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.008058 0.342255 0.564286 0.577749 0.838007 0.975213

# All cases

#pval_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) x[, 4])
#pval_icd_10_sh_i <- which(pval_icd_10_sh < thr) # J459  K573  K219  I251  I259  K529  M171  K297  M545  M179   K20  K449  I209 M199  F329   I10  M758  M159  I258  E119  J449 M4782  E039 M1999  E669  E780 F171  B968
#summary(unlist(pval_icd_10_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.000000 0.002039 0.040481 0.176427 0.265518 0.998576
#icd_vs_10_sh_sign <- do.call(rbind.data.frame, if_10_sh_vs_icd[pval_icd_10_sh_i])
#readme_icd_vs_10_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, not-combined, chapters 1-17, as dependent variable from binary coded SH PRS (10th percentile - 1, other - 0), whole sample of nonrelatives, cases"
#save(icd_vs_10_sh_sign, readme_icd_vs_10_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_ch1-17_vs_10_sh_all_nonr_cases_signif.RData")


# Test sample

#pval_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) x[, 4])
#pval_icd_10_sh_i <- which(pval_icd_10_sh < thr) # I10
#summary(unlist(pval_icd_10_sh))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0000013 0.0579529 0.2621130 0.3374100 0.5573738 0.9877788
#icd_vs_10_sh_sign <- do.call(rbind.data.frame, if_10_sh_vs_icd[pval_icd_10_sh_i])
#readme_icd_vs_10_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, not-combined, chapters 1-17, as dependent variable from binary coded SH PRS (10th percentile - 1, other - 0), test sample of nonrelatives, cases"
#save(icd_vs_10_sh_sign, readme_icd_vs_10_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_ch1-17_vs_10_sh_test_nonr_cases_signif.RData")

or_icd_10_10 <- lapply(if_10_10_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_10_10))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000001 0.499784 0.755625 0.797829 1.036627 3.394748
beta_icd_10_10 <- lapply(if_10_10_vs_icd, function(x) abs(x[, 1]))
or_icd_10_10_i <- which(beta_icd_10_10 > log(2))

# All cases

#or_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3144  0.6725  0.7578  0.7811  0.8722  1.5900
#beta_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) abs(x[, 1]))
#or_icd_10_sh_i <- which(beta_icd_10_sh > log(2)) 
#intersect(pval_icd_10_sh_i, or_icd_10_sh_i) # M758 M4782


# Test sample

#or_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3203  0.7025  0.8148  0.8339  0.9326  1.5842
#beta_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) abs(x[, 1]))
#or_icd_10_sh_i <- which(beta_icd_10_sh > log(2)) # J348  K29 I517
#intersect(pval_icd_10_sh_i, or_icd_10_sh_i) # 0 elements

# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.1 quantile for SH and 0.9 quantile for BP-SH PRS
pval_icd_10_90 <- lapply(if_10_90_vs_icd, function(x) x[, 4])
pval_icd_10_90_i <- which(pval_icd_10_90 < thr) # 0 elements
summary(unlist(pval_icd_10_90))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.01525 0.43433 0.67483 0.64219 0.90381 0.99909

# All cases

#pval_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) x[, 4])
#pval_icd_10_bp_sh_i <- which(pval_icd_10_bp_sh < thr) # J459 M179 N394 I209  I10 M139
#summary(unlist(pval_icd_10_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00000 0.03218 0.20962 0.33371 0.62932 0.99853
#icd_vs_10_bp_sh_sign <- do.call(rbind.data.frame, if_10_bp_sh_vs_icd[pval_icd_10_bp_sh_i])
#readme_icd_vs_10_bp_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, not-combined, chapters 1-17, as dependent variable from binary coded BP-SH PRS (10th percentile - 1, other - 0), whole sample of nonrelatives, cases"
#save(icd_vs_10_bp_sh_sign, readme_icd_vs_10_bp_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_ch1-17_vs_10_bp_sh_all_nonr_cases_signif.RData")


# Test sample

#pval_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) x[, 4])
#pval_icd_10_bp_sh_i <- which(pval_icd_10_bp_sh < thr) # 0 elements
#summary(unlist(pval_icd_10_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0005579 0.2160958 0.5131061 0.4774730 0.7023528 0.9844495

or_icd_10_90 <- lapply(if_10_90_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_10_90))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000013 0.5311088 0.7873624 0.8022280 1.0811379 2.8043476
beta_icd_10_90 <- lapply(if_10_90_vs_icd, function(x) abs(x[, 1]))
or_icd_10_90_i <- which(beta_icd_10_90 > log(2))

# All cases

#or_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7205  1.0148  1.1302  1.1351  1.2314  1.7836
#beta_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) abs(x[, 1]))
#or_icd_10_bp_sh_i <- which(beta_icd_10_bp_sh > log(2))
#intersect(pval_icd_10_bp_sh_i, or_icd_10_bp_sh_i) # 0

# Test sample

#or_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3525  0.8676  1.0047  1.0017  1.1099  1.8075
#beta_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) abs(x[, 1]))
#or_icd_10_bp_sh_i <- which(beta_icd_10_bp_sh > log(2))


# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.9 quantile for SH and 0.1 quantile for BP-SH PRS
pval_icd_90_10 <- lapply(if_90_10_vs_icd, function(x) x[, 4])
pval_icd_90_10_i <- which(pval_icd_90_10 < thr) # 0 elements
summary(unlist(pval_icd_90_10))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   0.0008202 0.3173015 0.5663905 0.5561854 0.8240230 0.9901085

or_icd_90_10 <- lapply(if_90_10_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_90_10))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000001 0.684244 1.063892 1.116819 1.459910 3.574734
beta_icd_90_10 <- lapply(if_90_10_vs_icd, function(x) abs(x[, 1]))
or_icd_90_10_i <- which(beta_icd_90_10 > log(2))

# All sample

#pval_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) x[, 4])
#pval_icd_90_sh_i <- which(pval_icd_90_sh < thr) # I849  M549  J459  K429  K573  M161  K219  M169  K259  N390 M2556  M754   K30 K802  I251  L031  K801 M2322  I259  K529  K625  M171  K297  G560  M545 M0699 M179   K20  K449  M751  J342  I209 M7966  M513  J440   J22  K590  K579  K519 G473  K589  K269  M796 M1991  G439  K760  M199  F329   I10 M5422  K921  J441 M159  K298  I258  E668  I739  I252  E119  J449 M4782 M2555  M478  M542 M1399 M5499 M5457  M069  E039  G409  M139  J439 M1999 M4796  M479 M1390  E669 M4799 E780  F171  E785
#summary(unlist(pval_icd_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000000 0.0000005 0.0029983 0.1183348 0.1128083 0.9751982

#icd_vs_90_sh_sign <- do.call(rbind.data.frame, if_90_sh_vs_icd[pval_icd_90_sh_i])
#readme_icd_vs_90_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, not-combined, chapters 1-17, as dependent variable from binary coded SH PRS (90th percentile - 1, other - 0), whole sample of nonrelatives, cases"
#save(icd_vs_90_sh_sign, readme_icd_vs_90_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_ch1-17_vs_90_sh_all_nonr_cases_signif.RData")

#or_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.5864  1.2128  1.3970  1.4116  1.5765  2.4002
#beta_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) abs(x[, 1]))
#or_icd_90_sh_i <- which(beta_icd_90_sh > log(2))
#intersect(pval_icd_90_sh_i, or_icd_90_sh_i) # M1991 K760 E668 M4782 M1390 M4799


# Test sample

#pval_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) x[, 4])
#pval_icd_90_sh_i <- which(pval_icd_90_sh < thr) # K219 M754 M545 M179 K449 F171
#summary(unlist(pval_icd_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000004 0.0439470 0.3050774 0.3417980 0.5706464 0.9807432

#icd_vs_90_sh_sign <- do.call(rbind.data.frame, if_90_sh_vs_icd[pval_icd_90_sh_i])
#readme_icd_vs_90_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, not-combined, chapters 1-17, as dependent variable from binary coded SH PRS (90th percentile - 1, other - 0), test sample of nonrelatives, cases"
#save(icd_vs_90_sh_sign, readme_icd_vs_90_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_ch1-17_vs_90_sh_test_nonr_cases_signif.RData")

#or_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.652   1.009   1.167   1.181   1.351   1.919
#beta_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) abs(x[, 1]))
#or_icd_90_sh_i <- which(beta_icd_90_sh > log(2))
#intersect(pval_icd_90_sh_i, or_icd_90_sh_i) # 0


# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.9 quantile for both SH and BP-SH PRS
pval_icd_90_90 <- lapply(if_90_90_vs_icd, function(x) x[, 4])
pval_icd_90_90_i <- which(pval_icd_90_90 < thr) # 0 elements
summary(unlist(pval_icd_90_90))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.000548 0.290169 0.544479 0.535535 0.810072 0.999354

or_icd_90_90 <- lapply(if_90_90_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_90_90))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000001 0.615884 1.016370 1.137191 1.510756 5.094163
beta_icd_90_90 <- lapply(if_90_90_vs_icd, function(x) abs(x[, 1]))
or_icd_90_90_i <- which(beta_icd_90_90 > log(2))

# All cases

#pval_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) x[, 4])
#pval_icd_90_bp_sh_i <- which(pval_icd_90_bp_sh < thr) # M171  M179  M199 M4782  M139 M1999
#summary(unlist(pval_icd_90_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.03192 0.25040 0.33095 0.57367 0.99507

#icd_vs_90_bp_sh_sign <- do.call(rbind.data.frame, if_90_bp_sh_vs_icd[pval_icd_90_bp_sh_i])
#readme_icd_vs_90_bp_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, not-combined, chapters 1-17, as dependent variable from binary coded BP-SH PRS (90th percentile - 1, other - 0), whole sample of nonrelatives, cases"
#save(icd_vs_90_bp_sh_sign, readme_icd_vs_90_bp_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_ch1-17_vs_90_bp_sh_all_nonr_cases_signif.RData")

#or_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.4233  0.7726  0.8690  0.8645  0.9565  1.3699
#beta_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) abs(x[, 1]))
#or_icd_90_bp_sh_i <- which(beta_icd_90_bp_sh > log(2))
#intersect(pval_icd_90_bp_sh_i, or_icd_90_bp_sh_i) # M4782


# Test sample

#pval_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) x[, 4])
#pval_icd_90_bp_sh_i <- which(pval_icd_90_bp_sh < thr) # 0
#summary(unlist(pval_icd_90_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.004072 0.284498 0.508252 0.521366 0.776025 0.990672

#or_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4502  0.8938  1.0031  1.0015  1.1128  1.7410
#beta_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) abs(x[, 1]))
#or_icd_90_bp_sh_i <- which(beta_icd_90_bp_sh > log(2))


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.1 quantile for both SH and BP-SH PRS
pval_opcs_10_10 <- lapply(if_10_10_vs_opcs, function(x) x[, 4])
pval_opcs_10_10_i <- which(pval_opcs_10_10 < thr) # 0 elements
summary(unlist(pval_opcs_10_10))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.02282 0.32957 0.61710 0.58499 0.83316 0.99432

or_opcs_10_10 <- lapply(if_10_10_vs_opcs, function(x) exp(x[, 1]))
summary(unlist(or_opcs_10_10))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  0.4798  0.7266  0.7966  1.1453  2.6146
beta_opcs_10_10 <- lapply(if_10_10_vs_opcs, function(x) abs(x[, 1]))
or_opcs_10_10_i <- which(beta_opcs_10_10 > log(2))

# All cases

#pval_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) x[, 4])
#pval_opcs_10_sh_i <- which(pval_opcs_10_sh < thr) # M459 G459 H259 G451 W822 H221 V544 A651 W903 W401 K633
#summary(unlist(pval_opcs_10_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.01003 0.07230 0.21894 0.38353 0.96079

#opcs_vs_10_sh_sign <- do.call(rbind.data.frame, if_10_sh_vs_opcs[pval_opcs_10_sh_i])
#readme_opcs_vs_10_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, not-combined, chapters X, Y, Z excluded, as dependent variable from binary coded SH PRS (10th percentile - 1, other - 0), whole sample of nonrelatives, cases"
#save(opcs_vs_10_sh_sign, readme_opcs_vs_10_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_10_sh_all_nonr_cases_signif.RData")

#or_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4204  0.7051  0.8125  0.8188  0.9109  1.5174
#beta_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) abs(x[, 1]))
#or_opcs_10_sh_i <- which(beta_opcs_10_sh > log(2))
#intersect(or_opcs_10_sh_i, pval_opcs_10_sh_i) # 0


# Test sample

#pval_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) x[, 4])
#pval_opcs_10_sh_i <- which(pval_opcs_10_sh < thr) # G451
#summary(unlist(pval_opcs_10_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000164 0.1272975 0.3662814 0.3917125 0.5751957 0.9749227

#opcs_vs_10_sh_sign <- do.call(rbind.data.frame, if_10_sh_vs_opcs[pval_opcs_10_sh_i])
#readme_opcs_vs_10_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, not-combined, chapters X, Y, Z excluded, as dependent variable from binary coded SH PRS (10th percentile - 1, other - 0), test sample of nonrelatives, cases"
#save(opcs_vs_10_sh_sign, readme_opcs_vs_10_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_10_sh_test_nonr_cases_signif.RData")

#or_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4816  0.7676  0.8585  0.8871  0.9469  2.9243
#beta_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) abs(x[, 1]))
#or_opcs_10_sh_i <- which(beta_opcs_10_sh > log(2))
#intersect(or_opcs_10_sh_i, pval_opcs_10_sh_i) # 0


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.1 quantile for SH and 0.9 uantile for BP-SH PRS
pval_opcs_10_90 <- lapply(if_10_90_vs_opcs, function(x) x[, 4])
pval_opcs_10_90_i <- which(pval_opcs_10_90 < thr) # 0 elements
summary(unlist(pval_opcs_10_90))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.1093  0.4409  0.6948  0.6622  0.9275  0.9961

or_opcs_10_90 <- lapply(if_10_90_vs_opcs, function(x) exp(x[, 1]))
summary(unlist(or_opcs_10_90))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  0.5327  0.8516  0.8521  1.1811  2.3389
beta_opcs_10_90 <- lapply(if_10_90_vs_opcs, function(x) abs(x[, 1]))
or_opcs_10_90_i <- which(beta_opcs_10_90 > log(2))

# All cases

#pval_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) x[, 4])
#pval_opcs_10_bp_sh_i <- which(pval_opcs_10_bp_sh < thr) # G451 W822 W401
#summary(unlist(pval_opcs_10_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.09711 0.36409 0.39232 0.66019 0.98754

#opcs_vs_10_bp_sh_sign <- do.call(rbind.data.frame, if_10_bp_sh_vs_opcs[pval_opcs_10_bp_sh_i])
#readme_opcs_vs_10_bp_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, not-combined, chapters X, Y, Z excluded, as dependent variable from binary coded BP-SH PRS (10th percentile - 1, other - 0), whole sample of nonrelatives, cases"
#save(opcs_vs_10_bp_sh_sign, readme_opcs_vs_10_bp_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_10_bp_sh_all_nonr_cases_signif.RData")

#or_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7169  0.9627  1.0730  1.0874  1.1731  1.5909
#beta_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) abs(x[, 1]))
#or_opcs_10_bp_sh_i <- which(beta_opcs_10_bp_sh > log(2))
#intersect(pval_opcs_10_bp_sh_i,or_opcs_10_bp_sh_i) # 0


# Test sample

#pval_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) x[, 4])
#pval_opcs_10_bp_sh_i <- which(pval_opcs_10_bp_sh < thr) # 0 elements
#summary(unlist(pval_opcs_10_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00898 0.17382 0.47021 0.46487 0.75514 0.97118

#or_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3244  0.8230  0.9464  0.9325  1.0487  1.5693
#beta_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) abs(x[, 1]))
#or_opcs_10_bp_sh_i <- which(beta_opcs_10_bp_sh > log(2))


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.9 quantile for SH and 0.1 uantile for BP-SH PRS
pval_opcs_90_10 <- lapply(if_90_10_vs_opcs, function(x) x[, 4])
pval_opcs_90_10_i <- which(pval_opcs_90_10 < thr) # 1 element; E25 3.601985e-05
summary(unlist(pval_opcs_90_10))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.000036 0.359297 0.598320 0.577755 0.841422 0.995399

or_opcs_90_10 <- lapply(if_90_10_vs_opcs, function(x) exp(x[, 1]))
summary(unlist(or_opcs_90_10))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  0.5327  0.8516  0.8521  1.1811  2.3389
beta_opcs_90_10 <- lapply(if_90_10_vs_opcs, function(x) abs(x[, 1]))
or_opcs_90_10_i <- which(beta_opcs_90_10 > log(2))

prs_90_10_opcs_i <- intersect(pval_opcs_90_10_i, or_opcs_90_10_i) # 1 element; E25 p = 3.601985e-05, or = 3.937817

# All cases

#pval_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) x[, 4])
#pval_opcs_90_sh_i <- which(pval_opcs_90_sh < thr) # H251 W852 E036 M459 G459 H259 G451 H229 W822 H221 V544 J183 E259 K634 W879 A651 T791 W903 A577 U212 A522 U051 W401 K633 W844 O291 S571 G199 U052 G211 A735 E369 V552 V559
#summary(unlist(pval_opcs_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000000 0.0000861 0.0385809 0.2034584 0.3536642 0.9860471

#opcs_vs_90_sh_sign <- do.call(rbind.data.frame, if_90_sh_vs_opcs[pval_opcs_90_sh_i])
#readme_opcs_vs_90_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, not-combined, chapters X, Y, Z excluded, as dependent variable from binary coded SH PRS (90th percentile - 1, other - 0), whole sample of nonrelatives, cases"
#save(opcs_vs_90_sh_sign, readme_opcs_vs_90_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_90_sh_all_nonr_cases_signif.RData")

#or_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.6231  1.1121  1.2778  1.3016  1.4669  2.4280
#beta_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) abs(x[, 1]))
#or_opcs_90_sh_i <- which(beta_opcs_90_sh > log(2))

# prs_90_sh_opcs_i <- intersect(pval_opcs_90_sh_i, or_opcs_90_sh_i) # E259 E369


# Test sample

#pval_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) x[, 4])
#pval_opcs_90_sh_i <- which(pval_opcs_90_sh < thr) # G459 G451 H229 W401 
#summary(unlist(pval_opcs_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000  0.1062  0.3175  0.3988  0.6731  0.9880

#opcs_vs_90_sh_sign <- do.call(rbind.data.frame, if_90_sh_vs_opcs[pval_opcs_90_sh_i])
#readme_opcs_vs_90_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, not-combined, chapters X, Y, Z excluded, as dependent variable from binary coded SH PRS (90th percentile - 1, other - 0), test sample of nonrelatives, cases"
#save(opcs_vs_90_sh_sign, readme_opcs_vs_90_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_90_sh_test_nonr_cases_signif.RData")

#or_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.3534  0.9920  1.1162  1.1436  1.3368  1.9773
#beta_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) abs(x[, 1]))
#or_opcs_90_sh_i <- which(beta_opcs_90_sh > log(2))

# prs_90_sh_opcs_i <- intersect(pval_opcs_90_sh_i, or_opcs_90_sh_i) # 0


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.9 quantile both for SH and BP-SH PRS
pval_opcs_90_90 <- lapply(if_90_90_vs_opcs, function(x) x[, 4])
pval_opcs_90_90_i <- which(pval_opcs_90_90 < thr) # 0 elements
summary(unlist(pval_opcs_90_90))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#  0.003216 0.338258 0.662570 0.601369 0.868501 0.996688

or_opcs_90_90 <- lapply(if_90_90_vs_opcs, function(x) exp(x[, 1]))
summary(unlist(or_opcs_90_90))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  0.6585  0.9715  1.0300  1.3699  3.5932
beta_opcs_90_90 <- lapply(if_90_90_vs_opcs, function(x) abs(x[, 1]))
or_opcs_90_90_i <- which(beta_opcs_90_90 > log(2))

# All cases

#pval_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) x[, 4])
#pval_opcs_90_bp_sh_i <- which(pval_opcs_90_bp_sh < thr) # G459 G451 W822 W401
#summary(unlist(pval_opcs_90_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.08858 0.32686 0.39841 0.67597 0.99746

#opcs_vs_90_bp_sh_sign <- do.call(rbind.data.frame, if_90_bp_sh_vs_opcs[pval_opcs_90_bp_sh_i])
#readme_opcs_vs_90_bp_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, not-combined, chapters X, Y, Z excluded, as dependent variable from binary coded BP-SH PRS (90th percentile - 1, other - 0), whole sample of nonrelatives, cases"
#save(opcs_vs_90_bp_sh_sign, readme_opcs_vs_90_bp_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_90_bp_sh_all_nonr_cases_signif.RData")

#or_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4736  0.8177  0.9236  0.9183  0.9904  1.5165
#beta_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) abs(x[, 1]))
#or_opcs_90_bp_sh_i <- which(beta_opcs_90_bp_sh > log(2))
#intersect(pval_opcs_90_bp_sh_i, or_opcs_90_bp_sh_i) # 0

# Test sample

#pval_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) x[, 4])
#pval_opcs_90_bp_sh_i <- which(pval_opcs_90_bp_sh < thr) # 0 elements
#summary(unlist(pval_opcs_90_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.004296 0.265883 0.537836 0.527552 0.798690 0.998465

#or_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2729  0.8943  0.9929  0.9979  1.1041  1.6843
#beta_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) abs(x[, 1]))
#or_opcs_90_bp_sh_i <- which(beta_opcs_90_bp_sh > log(2))

