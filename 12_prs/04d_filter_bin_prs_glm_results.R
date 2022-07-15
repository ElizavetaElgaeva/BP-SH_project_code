# Aim of this script is to filter the glm() results 
# for binary coded PRS traits against ICD10 and OPCS4 level 2 codes
# Note: test sample, non-relatieves only; ICD10 codes from chapters I-XVII only; no X, Y, Z OPCS codes

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load ICD10 glm results
load("icd10_level_2_chapter_1-17_vs_bin_prs_10_10.RData")
load("icd10_level_2_chapter_1-17_vs_bin_prs_10_90.RData")
load("icd10_level_2_chapter_1-17_vs_bin_sh_prs_90_10.RData")
load("icd10_level_2_chapter_1-17_vs_bin_prs_90_90.RData")

#load("/mnt/polyomica/projects/bp-sh/data/12_prs/icd10_level_2_chapter_1-17_vs_bin_sh_prs_90_sh.RData")
 
# Load OPCS glm results
load("opcs_level_2_no_xyz_vs_bin_prs_10_10.RData")
load("opcs_level_2_no_xyz_vs_bin_prs_10_90.RData")
load("opcs_level_2_no_xyz_vs_bin_sh_prs_90_10.RData")
load("opcs_level_2_no_xyz_vs_bin_prs_90_90.RData")

ls()

# Set significance threshold
thr <- 0.05 / (4*(length(if_10_10_vs_icd) + length(if_10_10_vs_opcs))) # 4.208754e-05; length(if_10_10_vs_icd) = 165, length(if_10_10_vs_opcs) = 132

# Start with ICD10 against binary coded PRS trait reflecting relatedness to 0.1 quantile for both SH and BP-SH PRS 
pval_icd_10_10 <- lapply(if_10_10_vs_icd, function(x) x[, 4])
pval_icd_10_10_i <- which(pval_icd_10_10 < thr) # 0 elements
summary(unlist(pval_icd_10_10))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.008058 0.342255 0.564286 0.577749 0.838007 0.975213

#pval_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) x[, 4])
#pval_icd_10_sh_i <- which(pval_icd_10_sh < thr) # I10, J45, K29
#summary(unlist(pval_icd_10_sh))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0000013 0.0592227 0.2223695 0.3054859 0.4854355 0.9877788

or_icd_10_10 <- lapply(if_10_10_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_10_10))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000001 0.499784 0.755625 0.797829 1.036627 3.394748
beta_icd_10_10 <- lapply(if_10_10_vs_icd, function(x) abs(x[, 1]))
or_icd_10_10_i <- which(beta_icd_10_10 > log(2))

#or_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3653  0.6988  0.7918  0.8303  0.9303  1.7909
#beta_icd_10_sh <- lapply(if_10_sh_vs_icd, function(x) abs(x[, 1]))
#or_icd_10_sh_i <- which(beta_icd_10_sh > log(2)) # E05, I77, J34, J38, K82, L57
#intersect(pval_icd_10_sh, or_icd_10_sh) # 0 elements

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
#0.006336 0.229814 0.495463 0.495172 0.732901 0.998132

or_icd_10_90 <- lapply(if_10_90_vs_icd, function(x) exp(x[, 1]))
summary(unlist(or_icd_10_90))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000013 0.5311088 0.7873624 0.8022280 1.0811379 2.8043476
beta_icd_10_90 <- lapply(if_10_90_vs_icd, function(x) abs(x[, 1]))
or_icd_10_90_i <- which(beta_icd_10_90 > log(2))


#or_icd_10_bp_sh <- lapply(if_10_bp_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4414  0.8848  0.9874  1.0000  1.1123  1.5940
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
#pval_icd_90_sh_i <- which(pval_icd_90_sh < thr) # F17, K21, K29, K44, M13, M17, M19, M54
#summary(unlist(pval_icd_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.00000 0.02984 0.25011 0.34865 0.63762 0.99808

#or_icd_90_sh <- lapply(if_90_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.5325  1.0069  1.1821  1.1951  1.3636  2.1124
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
#0.00421 0.23208 0.52540 0.50816 0.75195 0.98618

#or_icd_90_bp_sh <- lapply(if_90_bp_sh_vs_icd, function(x) exp(x[, 1]))
#summary(unlist(or_icd_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.4350  0.8857  0.9745  0.9885  1.0896  1.5418
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
#pval_opcs_10_sh_i <- which(pval_opcs_10_sh < thr) # G45
#summary(unlist(pval_opcs_10_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.0000002 0.1112621 0.3416465 0.3921942 0.6297752 0.9792937

#or_opcs_10_sh <- lapply(if_10_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_10_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2435  0.7141  0.8582  0.8438  0.9559  1.6042
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
#0.009091 0.160707 0.379945 0.452032 0.754535 0.996036

#or_opcs_10_bp_sh <- lapply(if_10_bp_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_10_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3481  0.8054  0.9711  0.9691  1.1364  1.5795
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
#pval_opcs_90_sh_i <- which(pval_opcs_90_sh < thr) # A52, G45, T25
#summary(unlist(pval_opcs_90_sh))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.0000  0.1323  0.3888  0.4216  0.6852  0.9944

#or_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_90_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.5804  0.9654  1.1215  1.1604  1.2996  2.7033
#beta_opcs_90_sh <- lapply(if_90_sh_vs_opcs, function(x) abs(x[, 1]))
#or_opcs_90_sh_i <- which(beta_opcs_90_sh > log(2))

# prs_90_sh_opcs_i <- intersect(pval_opcs_90_sh_i, or_opcs_90_sh_i) # T25


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
#0.008549 0.278266 0.509688 0.508647 0.749071 0.996145

#or_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) exp(x[, 1]))
#summary(unlist(or_opcs_90_bp_sh))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4462  0.8926  1.0185  1.0328  1.1725  1.7489
#beta_opcs_90_bp_sh <- lapply(if_90_bp_sh_vs_opcs, function(x) abs(x[, 1]))
#or_opcs_90_bp_sh_i <- which(beta_opcs_90_bp_sh > log(2))

