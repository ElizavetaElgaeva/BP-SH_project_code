# Aim of this script is to filter the glm() results 
# for binary coded PRS traits against ICD10 and OPCS4 level 2 codes
# Note: test sample, non-relatieves only; ICD10 codes from chapters I-XVII only; no X, Y, Z OPCS codes

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load ICD10 glm results
load("icd10_level_2_chapter_1-17_vs_bin_prs_10_10.RData")
load("icd10_level_2_chapter_1-17_vs_bin_prs_10_90.RData")
load("icd10_level_2_chapter_1-17_vs_bin_sh_prs_90_10.RData")
load("icd10_level_2_chapter_1-17_vs_bin_prs_90_90.RData")
# load("icd10_vs_bin_sh_prs_90_sh.RData")
#load("/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter_1-17_vs_bin_sh_prs_90_sh.RData")
 
# Load OPCS glm results
load("opcs_level_2_no_xyz_vs_bin_prs_10_10.RData")
load("opcs_level_2_no_xyz_vs_bin_prs_10_90.RData")
load("opcs_level_2_no_xyz_vs_bin_sh_prs_90_10.RData")
load("opcs_level_2_no_xyz_vs_bin_prs_90_90.RData")

# load("opcs_vs_bin_prs_10_sh.RData")
#load("opcs_bin_prs_10_bp_sh.RData")
ls()

# Set significance threshold
thr <- 0.05 / (4*(length(if_10_10_vs_icd) + length(if_10_10_vs_opcs))) # 4.208754e-05; length(if_10_10_vs_icd) = 165, length(if_10_10_vs_opcs) = 132
#thr <- 0.05 / (4*(length(if_10_sh_vs_icd) + length(if_90_sh_vs_opcs))) # 2.73523e-05; length(if_10_sh_vs_icd) = 267, length(if_90_bp_sh_vs_opcs) = 190


# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.1 quantile for both SH and BP-SH PRS 
pval_icd_10_10 <- lapply(if_10_10_vs_icd, function(x) x[, 4])
pval_icd_10_10_i <- which(pval_icd_10_10 < thr) # 0 elements
summary(unlist(pval_icd_10_10))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.008058 0.342255 0.564286 0.577749 0.838007 0.975213

#pval_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) x[, 4])
#pval_icd_10_sh_i <- which(pval_icd_10_sh < thr) # I10
#summary(unlist(pval_icd_10_sh))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0000013 0.0595199 0.2203845 0.3285945 0.5591212 0.9877788

or_icd_10_10 <- lapply(if_10_10_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_10_10))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000001 0.499784 0.755625 0.797829 1.036627 3.394748
beta_icd_10_10 <- lapply(if_10_10_vs_icd, function(x) abs(x[, 1]))
or_icd_10_10_i <- which(beta_icd_10_10 > log(2))

#or_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3203  0.7121  0.8189  0.8516  0.9511  1.5842
#beta_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) abs(x[, 1]))
#or_icd_10_sh_i <- which(beta_icd_10_sh > log(2)) # J348  K29 I517
#intersect(pval_icd_10_sh_i, or_icd_10_sh_i) # 0 elements

# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.1 quantile for SH and 0.9 quantile for BP-SH PRS
pval_icd_10_90 <- lapply(if_10_90_vs_icd, function(x) x[, 4])
pval_icd_10_90_i <- which(pval_icd_10_90 < thr) # 0 elements
summary(unlist(pval_icd_10_90))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.01525 0.43433 0.67483 0.64219 0.90381 0.99909


#pval_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) x[, 4])
#pval_icd_10_bp_sh_i <- which(pval_icd_10_bp_sh < thr) # 0 elements
#summary(unlist(pval_icd_10_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0005579 0.2143139 0.5055283 0.4837844 0.7142593 0.9963940

or_icd_10_90 <- lapply(if_10_90_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_10_90))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000013 0.5311088 0.7873624 0.8022280 1.0811379 2.8043476
beta_icd_10_90 <- lapply(if_10_90_vs_icd, function(x) abs(x[, 1]))
or_icd_10_90_i <- which(beta_icd_10_90 > log(2))


#or_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2726  0.8777  1.0032  1.0011  1.1226  1.8075
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


#pval_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) x[, 4])
#pval_icd_90_sh_i <- which(pval_icd_90_sh < thr) # K219 M754 M545 M179 K449 F171
#summary(unlist(pval_icd_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000004 0.0504154 0.3091381 0.3568354 0.6085242 0.9882826

#or_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.652   1.006   1.156   1.165   1.321   1.919
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


#pval_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) x[, 4])
#pval_icd_90_bp_sh_i <- which(pval_icd_90_bp_sh < thr) # 0 elements
#summary(unlist(pval_icd_90_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.004072 0.296856 0.515964 0.520643 0.771527 0.990672

#or_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4502  0.8900  0.9980  1.0037  1.1201  1.7410
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


#pval_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) x[, 4])
#pval_opcs_10_sh_i <- which(pval_opcs_10_sh < thr) # X998 G451 Z943
#summary(unlist(pval_opcs_10_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000046 0.0937551 0.3626877 0.3785776 0.5853029 0.9979890

#or_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4288  0.7651  0.8542  0.8907  0.9534  2.9243
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


#pval_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) x[, 4])
#pval_opcs_10_bp_sh_i <- which(pval_opcs_10_bp_sh < thr) # 0 elements
#summary(unlist(pval_opcs_10_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00898 0.23841 0.46454 0.47583 0.74156 0.97118

#or_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3244  0.8393  0.9497  0.9507  1.0661  1.5693
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


#pval_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) x[, 4])
#pval_opcs_90_sh_i <- which(pval_opcs_90_sh < thr) # X998 G459 G451 H229 Z942 Z274 Y767 Z812
#summary(unlist(pval_opcs_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.08527 0.32319 0.39537 0.66525 0.98799

#or_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3534  0.9877  1.1162  1.1375  1.3119  2.0877
#beta_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) abs(x[, 1]))
#or_opcs_90_sh_i <- which(beta_opcs_90_sh > log(2))

# prs_90_sh_opcs_i <- intersect(pval_opcs_90_sh_i, or_opcs_90_sh_i) # Z812


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


#pval_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) x[, 4])
#pval_opcs_90_bp_sh_i <- which(pval_opcs_90_bp_sh < thr) # 0 elements
#summary(unlist(pval_opcs_90_bp_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.004296 0.254544 0.523665 0.514358 0.774792 0.998465

#or_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2729  0.9014  1.0009  1.0098  1.1160  1.6843
#beta_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) abs(x[, 1]))
#or_opcs_90_bp_sh_i <- which(beta_opcs_90_bp_sh > log(2))

