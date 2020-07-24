# Aim of this script is to estimate heritability of Back pain after
# shared heredity subtraction and its genetic correlation with shared heredity
# using LD Score implementd in GWAS-Map (from test container). Discovery study

# Estimate heritability for BP-SH
run_ldscore --h2 --gwas-id=41
# 0.0141 ## se 0.002 

# Estimate genetic correlation for SH and BP-SH
run_ldscore --rg --gwas-id=19,41
## -0.0031 (se 0.051)

# compare with the results of 01a script
