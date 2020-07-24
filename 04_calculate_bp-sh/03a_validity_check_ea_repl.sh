# Aim of this script is to estimate heritability of Back pain after
# shared heredity subtraction and its genetic correlation with shared heredity
# using LD Score implementd in GWAS-Map (from test container). EA Replication study

# Estimate heritability for BP-SH
run_ldscore --h2 --gwas-id=42
# 0.011 ## 0.003

# Estimate genetic correlation for SH and BP-SH
run_ldscore --rg --gwas-id=26,42
## -0.0117 (se 0.078)

# Estimate genetic correlation for BP-SH Discovery and EA Replication 
## 0.8567

# compare with the results of 01a script
