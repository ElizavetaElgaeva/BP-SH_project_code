# Aim of this script is to estimate heritability of
# shared heredity and its genetic correlations with 6 pain traits (AA Replication study)
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for SH (AA Replication stuudy)
run_ldscore --h2 --gwas-id=40
## h2=-0.0039  se=0.0486

# Estimate pairwise genetic correlations for SH and 6 pain traits (AA Replication study)
run_ldscore --rg --gwas-id=34,40 # sh and Hip pain
## None
run_ldscore --rg --gwas-id=35,40 # sh and Back pain
## None
run_ldscore --rg --gwas-id=36,40 # sh and Neck pain
##  None
run_ldscore --rg --gwas-id=37,40 # sh and Knee pain
##  None
run_ldscore --rg --gwas-id=38,40 # sh and Head pain
##  None
run_ldscore --rg --gwas-id=39,40 # sh and Stomach pain
##  None

# Estimate pairwise genetic correlation for SH from Discovery, EA, AA and SA Replication studies
run_ldscore --rg --gwas-id=19,40 # aa sh and discovery sh
## ERROR
run_ldscore --rg --gwas-id=26,40 # aa sh and ea sh
## ERROR
run_ldscore --rg --gwas-id=33,40 # aa sh and sa sh
## None



# compare with the results of 07a script
