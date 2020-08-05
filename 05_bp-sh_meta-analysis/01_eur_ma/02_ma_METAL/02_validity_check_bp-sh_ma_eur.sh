# Aim of this script is to estimate heritability of
# European meta-analysis of Back pain without shared heredity
# and its genetic correlations with Shared heredity and Back pain without
# Shared heredity traits from different datasets
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for BP-SH MA EUR
run_ldscore --h2 --gwas-id=49
## h2=0.0121 se=0.0013

# Estimate pairwise genetic correlations for BP-SH MA EUR
run_ldscore --rg --gwas-id=49,41 # bp-sh ma eur and bp-sh discovery
## 0.9796 se=0.022
run_ldscore --rg --gwas-id=49,42 # bp-sh ma eur and bp-sh eur replication
## 0.9491 se=0.069
run_ldscore --rg --gwas-id=49,47 # bp-sh ma eur and sh ma eur
## 4.5791e-07 se=0.039
run_ldscore --rg --gwas-id=49,48 # bp-sh ma eur and sh ma repl
## -0.0278 se=0.048
run_ldscore --rg --gwas-id=49,19 # bp-sh ma eur and sh discovery
## 0.0324 se=0.045
run_ldscore --rg --gwas-id=49,26 # bp-sh ma eur and sh eur replication
## -0.0231 se=0.044

