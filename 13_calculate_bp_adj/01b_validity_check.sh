# Aim of this script is to estimate heritability of CBP
# adjusted for all other original pain traits
# using LD Score implementd in GWAS-Map (from test container). Discovery study

# Estimate heritability for adjusted BP
run_ldscore --h2 --gwas-id=53
# 0.007 ## se 0.0018 

# Estimate genetic correlation for SH and adjusted BP
run_ldscore --rg --gwas-id=19,53
## 0.1537 (se 0.073, p 0.0351)

# Estimate genetic correlation for BP and adjusted BP
run_ldscore --rg --gwas-id=14,53
## 0.4912 (se 0.065, p 6.0025e-14)

# Estimate genetic correlation for BP-SH and adjusted BP
run_ldscore --rg --gwas-id=41,53
## 0.7986 (se 0.042, p 9.0068e-80)

# Estimate genetic correlation for Hip pain and adjusted BP
run_ldscore --rg --gwas-id=13,53
## 0.0253 (se 0.098, p 0.7958)

# Estimate genetic correlation for Neck pain and adjusted BP
run_ldscore --rg --gwas-id=15,53
## 0.0146 (se 0.108, p 0.8926)

# Estimate genetic correlation for Knee pain and adjusted BP
run_ldscore --rg --gwas-id=16,53
## 0.0165 (se 0.09, p 0.8549)

# Estimate genetic correlation for Head pain and adjusted BP
run_ldscore --rg --gwas-id=17,53
## -0.0018 (se 0.083, p 0.9828)

# Estimate genetic correlation for Stomach pain and adjusted BP
run_ldscore --rg --gwas-id=18,53
## -0.0033 (se 0.117, p 0.9775)

# compare with the results of 01a script
