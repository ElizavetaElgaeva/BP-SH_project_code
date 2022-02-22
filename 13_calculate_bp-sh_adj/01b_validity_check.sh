# Aim of this script is to estimate heritability of BP-SH
# adjusted for head and knee pain and its genetic correlations with other traits
# using LD Score implementd in GWAS-Map (from test container). Discovery study

# Estimate heritability for adjusted BP-SH
run_ldscore --h2 --gwas-id=51
# 0.0131 ## se 0.0018 

# Estimate genetic correlation for SH and adjusted BP-SH
run_ldscore --rg --gwas-id=19,51
## 0.3316 (se 0.05, p 2.6298e-11)

# Estimate genetic correlation for BP and adjusted BP-SH
run_ldscore --rg --gwas-id=14,51
## 0.7081 (se 0.029, p 7.4314e-133)

# Estimate genetic correlation for BP-SH and adjusted BP-SH
run_ldscore --rg --gwas-id=41,51
## 0.9305 (se 0.01, p 0.0)

# compare with the results of 01a script
