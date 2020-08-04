# Aim of this script is to estimate heritability of
# European meta-analysis of shared heredity and its genetic correlations with
# Shared heredity from Discovery study and European Replication study
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for SH MA EUR
run_ldscore --h2 --gwas-id=47
## h2=0.0734  se=0.0025

# Estimate pairwise genetic correlations for SH MA EUR
run_ldscore --rg --gwas-id=47,19 # sh ma eur and sh discovery
## 0.9995 se=0.005
run_ldscore --rg --gwas-id=47,26 # sh ma eur and sh eur replication
## 0.9909 se=0.01

