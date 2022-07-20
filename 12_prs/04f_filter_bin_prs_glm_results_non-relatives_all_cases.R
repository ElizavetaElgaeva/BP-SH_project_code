# Aim of this script is to filter the glm() results
# for binary coded PRS traits against ICD10 and OPCS4 not-combined
# Note: whole sample, cases, non-relatieves only; ICD10 codes from chapters I-XVII only; no X, Y, Z OPCS codes

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load ICD10 glm results
load("icd10_chapter1-17_vs_bin_prs_10_sh_all_cases.RData")
load("icd10_chapter1-17vs_bin_sh_prs_90_sh_all_cases.RData")
load("icd10_chapter1-17_vs_bin_prs_10_bp_sh_all_cases.RData")
load("icd10_chapter1-17_vs_bin_prs_90_bp_sh_all_cases.RData")

# Load OPCS glm results
load("opcs_no_xyz_vs_bin_prs_10_sh_all_cases.RData")
load("opcs_no_xyz_bin_prs_10_bp_sh_all_cases.RData")
load("opcs_no_xyz_vs_bin_prs_90_sh_all_cases.RData")
load("opcs_no_xyz_vs_bin_prs_90_bp_sh_all_cases.RData")

ls()

# Set significance threshold
thr <- 0.05 / (4*(length(if_10_sh_vs_icd) + length(if_90_sh_vs_opcs))) # 3.140704e-05; length(if_10_sh_vs_icd) = 243, length(if_90_sh_vs_opcs) = 155 (all nonrel cases)

# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.1 quantile of SH PRS 
pval_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) x[, 4])
pval_icd_10_sh_i <- which(pval_icd_10_sh < thr) # J459  K573  K219  I251  I259  K529  M171  K297  M545  M179   K20  K449  I209 M199  F329   I10  M758  M159  I258  E119  J449 M4782  E039 M1999  E669  E780 F171  B968
summary(unlist(pval_icd_10_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.000000 0.002039 0.040481 0.176427 0.265518 0.998576
icd_vs_10_sh_sign <- do.call(rbind.data.frame, if_10_sh_vs_icd[pval_icd_10_sh_i])
readme_icd_vs_10_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, not-combined, chapters 1-17, as dependent variable from binary coded SH PRS (10th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(icd_vs_10_sh_sign, readme_icd_vs_10_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_ch1-17_vs_10_sh_all_nonr_cases_signif.RData")

or_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3144  0.6725  0.7578  0.7811  0.8722  1.5900
beta_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) abs(x[, 1]))
or_icd_10_sh_i <- which(beta_icd_10_sh > log(2)) 
intersect(pval_icd_10_sh_i, or_icd_10_sh_i) # M758 M4782


# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.1 quantile of BP-SH PRS
pval_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) x[, 4])
pval_icd_10_bp_sh_i <- which(pval_icd_10_bp_sh < thr) # J459 M179 N394 I209  I10 M139
summary(unlist(pval_icd_10_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00000 0.03218 0.20962 0.33371 0.62932 0.99853
icd_vs_10_bp_sh_sign <- do.call(rbind.data.frame, if_10_bp_sh_vs_icd[pval_icd_10_bp_sh_i])
readme_icd_vs_10_bp_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, not-combined, chapters 1-17, as dependent variable from binary coded BP-SH PRS (10th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(icd_vs_10_bp_sh_sign, readme_icd_vs_10_bp_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_ch1-17_vs_10_bp_sh_all_nonr_cases_signif.RData")

or_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7205  1.0148  1.1302  1.1351  1.2314  1.7836
beta_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) abs(x[, 1]))
or_icd_10_bp_sh_i <- which(beta_icd_10_bp_sh > log(2))
intersect(pval_icd_10_bp_sh_i, or_icd_10_bp_sh_i) # 0


# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.9 quantile of SH PRS
pval_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) x[, 4])
pval_icd_90_sh_i <- which(pval_icd_90_sh < thr) # I849  M549  J459  K429  K573  M161  K219  M169  K259  N390 M2556  M754   K30 K802  I251  L031  K801 M2322  I259  K529  K625  M171  K297  G560  M545 M0699 M179   K20  K449  M751  J342  I209 M7966  M513  J440   J22  K590  K579  K519 G473  K589  K269  M796 M1991  G439  K760  M199  F329   I10 M5422  K921  J441 M159  K298  I258  E668  I739  I252  E119  J449 M4782 M2555  M478  M542 M1399 M5499 M5457  M069  E039  G409  M139  J439 M1999 M4796  M479 M1390  E669 M4799 E780  F171  E785
summary(unlist(pval_icd_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000000 0.0000005 0.0029983 0.1183348 0.1128083 0.9751982
icd_vs_90_sh_sign <- do.call(rbind.data.frame, if_90_sh_vs_icd[pval_icd_90_sh_i])
readme_icd_vs_90_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, not-combined, chapters 1-17, as dependent variable from binary coded SH PRS (90th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(icd_vs_90_sh_sign, readme_icd_vs_90_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_ch1-17_vs_90_sh_all_nonr_cases_signif.RData")

or_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.5864  1.2128  1.3970  1.4116  1.5765  2.4002
beta_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) abs(x[, 1]))
or_icd_90_sh_i <- which(beta_icd_90_sh > log(2))
intersect(pval_icd_90_sh_i, or_icd_90_sh_i) # M1991 K760 E668 M4782 M1390 M4799


# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.9 quantile of BP-SH PRS
pval_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) x[, 4])
pval_icd_90_bp_sh_i <- which(pval_icd_90_bp_sh < thr) # M171  M179  M199 M4782  M139 M1999
summary(unlist(pval_icd_90_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.03192 0.25040 0.33095 0.57367 0.99507
icd_vs_90_bp_sh_sign <- do.call(rbind.data.frame, if_90_bp_sh_vs_icd[pval_icd_90_bp_sh_i])
readme_icd_vs_90_bp_sh_sign <- "Statistically significant results of glm conidering ICD10 codes, not-combined, chapters 1-17, as dependent variable from binary coded BP-SH PRS (90th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(icd_vs_90_bp_sh_sign, readme_icd_vs_90_bp_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/icd_ch1-17_vs_90_bp_sh_all_nonr_cases_signif.RData")

or_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.4233  0.7726  0.8690  0.8645  0.9565  1.3699
beta_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) abs(x[, 1]))
or_icd_90_bp_sh_i <- which(beta_icd_90_bp_sh > log(2))
intersect(pval_icd_90_bp_sh_i, or_icd_90_bp_sh_i) # M4782


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.1 quantile of SH PRS
pval_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) x[, 4])
pval_opcs_10_sh_i <- which(pval_opcs_10_sh < thr) # M459 G459 H259 G451 W822 H221 V544 A651 W903 W401 K633
summary(unlist(pval_opcs_10_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.01003 0.07230 0.21894 0.38353 0.96079
opcs_vs_10_sh_sign <- do.call(rbind.data.frame, if_10_sh_vs_opcs[pval_opcs_10_sh_i])
readme_opcs_vs_10_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, not-combined, chapters X, Y, Z excluded, as dependent variable from binary coded SH PRS (10th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(opcs_vs_10_sh_sign, readme_opcs_vs_10_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_10_sh_all_nonr_cases_signif.RData")

or_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) exp(x[, 1]))
summary(unlist(or_opcs_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4204  0.7051  0.8125  0.8188  0.9109  1.5174
beta_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) abs(x[, 1]))
or_opcs_10_sh_i <- which(beta_opcs_10_sh > log(2))
intersect(or_opcs_10_sh_i, pval_opcs_10_sh_i) # 0


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.1 quantile of BP-SH PRS
pval_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) x[, 4])
pval_opcs_10_bp_sh_i <- which(pval_opcs_10_bp_sh < thr) # G451 W822 W401
summary(unlist(pval_opcs_10_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.09711 0.36409 0.39232 0.66019 0.98754
opcs_vs_10_bp_sh_sign <- do.call(rbind.data.frame, if_10_bp_sh_vs_opcs[pval_opcs_10_bp_sh_i])
readme_opcs_vs_10_bp_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, not-combined, chapters X, Y, Z excluded, as dependent variable from binary coded BP-SH PRS (10th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(opcs_vs_10_bp_sh_sign, readme_opcs_vs_10_bp_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_10_bp_sh_all_nonr_cases_signif.RData")

or_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) exp(x[, 1]))
summary(unlist(or_opcs_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.7169  0.9627  1.0730  1.0874  1.1731  1.5909
beta_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) abs(x[, 1]))
or_opcs_10_bp_sh_i <- which(beta_opcs_10_bp_sh > log(2))
intersect(pval_opcs_10_bp_sh_i,or_opcs_10_bp_sh_i) # 0


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.9 quantile of SH PRS
pval_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) x[, 4])
pval_opcs_90_sh_i <- which(pval_opcs_90_sh < thr) # H251 W852 E036 M459 G459 H259 G451 H229 W822 H221 V544 J183 E259 K634 W879 A651 T791 W903 A577 U212 A522 U051 W401 K633 W844 O291 S571 G199 U052 G211 A735 E369 V552 V559
summary(unlist(pval_opcs_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000000 0.0000861 0.0385809 0.2034584 0.3536642 0.9860471
opcs_vs_90_sh_sign <- do.call(rbind.data.frame, if_90_sh_vs_opcs[pval_opcs_90_sh_i])
readme_opcs_vs_90_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, not-combined, chapters X, Y, Z excluded, as dependent variable from binary coded SH PRS (90th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(opcs_vs_90_sh_sign, readme_opcs_vs_90_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_90_sh_all_nonr_cases_signif.RData")

or_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) exp(x[, 1]))
summary(unlist(or_opcs_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.6231  1.1121  1.2778  1.3016  1.4669  2.4280
beta_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) abs(x[, 1]))
or_opcs_90_sh_i <- which(beta_opcs_90_sh > log(2))
intersect(pval_opcs_90_sh_i, or_opcs_90_sh_i) # E259 E369


# Start with OPCS against binary coded PRS trait reflecting relatedness to 0.9 quantile both of BP-SH PRS
pval_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) x[, 4])
pval_opcs_90_bp_sh_i <- which(pval_opcs_90_bp_sh < thr) # G459 G451 W822 W401
summary(unlist(pval_opcs_90_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.08858 0.32686 0.39841 0.67597 0.99746
opcs_vs_90_bp_sh_sign <- do.call(rbind.data.frame, if_90_bp_sh_vs_opcs[pval_opcs_90_bp_sh_i])
readme_opcs_vs_90_bp_sh_sign <- "Statistically significant results of glm conidering OPCS4 codes, not-combined, chapters X, Y, Z excluded, as dependent variable from binary coded BP-SH PRS (90th percentile - 1, other - 0), whole sample of nonrelatives, cases"
save(opcs_vs_90_bp_sh_sign, readme_opcs_vs_90_bp_sh_sign, file = "/mnt/polyomica/projects/bp-sh/data/12_prs/opcs_no_xyz_vs_90_bp_sh_all_nonr_cases_signif.RData")

or_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) exp(x[, 1]))
summary(unlist(or_opcs_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4736  0.8177  0.9236  0.9183  0.9904  1.5165
beta_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) abs(x[, 1]))
or_opcs_90_bp_sh_i <- which(beta_opcs_90_bp_sh > log(2))
intersect(pval_opcs_90_bp_sh_i, or_opcs_90_bp_sh_i) # 0


