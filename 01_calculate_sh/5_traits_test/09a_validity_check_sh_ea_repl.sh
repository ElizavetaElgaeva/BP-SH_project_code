# Aim of this script is to estimate heritability of
# shared heredity and its genetic correlations with 6 pain traits (EA Replication study)
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for SH (EA Replication stuudy)
run_ldscore --h2 --gwas-id=26
## h2=0.0796  se=0.0041

# Estimate pairwise genetic correlations for SH and 6 pain traits (EA Replication study)
run_ldscore --rg --gwas-id=20,26 # sh and Hip pain
## 0.9055 se=0.04
run_ldscore --rg --gwas-id=21,26 # sh and Back pain
##  0.9214 se=0.018
run_ldscore --rg --gwas-id=22,26 # sh and Neck pain
##  0.9097 se=0.019
run_ldscore --rg --gwas-id=23,26 # sh and Knee pain
##  0.8482 se=0.021
run_ldscore --rg --gwas-id=24,26 # sh and Head pain
##  0.545 se=0.036
run_ldscore --rg --gwas-id=25,26 # sh and Stomach pain
##  0.7673 se=0.043

# Estimate pairwise genetic correlation for SH from Discovery and EA Replication studies
run_ldscore --rg --gwas-id=19,26
##  0.9836 se=0.029

# compare with the results of 07a script
