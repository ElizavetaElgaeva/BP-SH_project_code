# Aim of this script is to estimate heritability of
# shared heredity and its genetic correlations with 6 pain traits (Discovery study)
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for SH (Discovery stuudy)
run_ldscore --h2 --gwas-id=19
## h2=0.0706  se=0.0034

# Estimate pairwise genetic correlations for SH and 6 pain traits (Discovery study)
run_ldscore --rg --gwas-id=13,19 # sh and Hip pain
## 0.881 se=0.025
run_ldscore --rg --gwas-id=14,19 # sh and Back pain
##  0.8991 se=0.012
run_ldscore --rg --gwas-id=15,19 # sh and Neck pain
##  0.9403 se=0.016
run_ldscore --rg --gwas-id=16,19 # sh and Knee pain
##  0.7826 se=0.019
run_ldscore --rg --gwas-id=17,19 # sh and Head pain
##  0.5399 se=0.026
run_ldscore --rg --gwas-id=18,19 # sh and Stomach pain
##  0.8311 se=0.04


# compare with the results of 07a script
