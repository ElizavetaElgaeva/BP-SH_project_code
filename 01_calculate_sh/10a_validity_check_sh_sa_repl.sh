# Aim of this script is to estimate heritability of
# shared heredity and its genetic correlations with 6 pain traits (SA Replication study)
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for SH (SA Replication stuudy)
run_ldscore --h2 --gwas-id=33
## h2=0.0289  se=0.0528

# Estimate pairwise genetic correlations for SH and 6 pain traits (SA Replication study)
run_ldscore --rg --gwas-id=27,33 # sh and Hip pain
## 0.9725 se=0.604
run_ldscore --rg --gwas-id=28,33 # sh and Back pain
## None
run_ldscore --rg --gwas-id=29,33 # sh and Neck pain
##  1.3631 se=0.747
run_ldscore --rg --gwas-id=30,33 # sh and Knee pain
##  0.751 se=0.522
run_ldscore --rg --gwas-id=31,33 # sh and Head pain
##  0.1371 se=0.695
run_ldscore --rg --gwas-id=32,33 # sh and Stomach pain
##  ERROR

# Estimate pairwise genetic correlation for SH from Discovery, EA and SA Replication studies
run_ldscore --rg --gwas-id=19,33 # sa sh and discovery sh
##  0.8871 se=0.722
run_ldscore --rg --gwas-id=26,33 # sa sh and ea sh
## 1.2559 se=0.973


# compare with the results of 07a script
