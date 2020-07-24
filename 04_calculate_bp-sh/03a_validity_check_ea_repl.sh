# Aim of this script is to estimate heritability of Back pain after
# shared heredity subtraction and its genetic correlation with shared heredity
# using LD Score implementd in GWAS-Map (from test container). EA Replication study

# Estimate heritability for BP-SH
run_ldscore --h2 --gwas-id=42
#  ## 

# Estimate genetic correlation for SH and BP-SH
run_ldscore --rg --gwas-id=26,42
##  (se )

# compare with the results of 01a script
