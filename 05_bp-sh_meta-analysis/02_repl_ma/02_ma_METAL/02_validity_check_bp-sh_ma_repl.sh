# Aim of this script is to estimate heritability of
# Replication meta-analysis of Back pain without shared heredity and its genetic correlations with
# Shared heredity and Back pain without shared heredity from different studies
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for BP-SH MA Repl
run_ldscore --h2 --gwas-id=50
## h2=0.0104 se=0.0024

# Estimate pairwise genetic correlations for BP-SH MA Repl
run_ldscore --rg --gwas-id=50,19 # bp-sh ma repl and sh discovery
## 0.0974 se=0.063
run_ldscore --rg --gwas-id=50,26 # bp-sh ma repl and sh eur replication
## -0.0106 se=0.078
run_ldscore --rg --gwas-id=50,47 # bp-sh ma repl and sh ma eur
## 0.0399 se=0.054
run_ldscore --rg --gwas-id=50,48 # bp-sh ma repl and sh ma repl
## -0.0204 se=0.072
run_ldscore --rg --gwas-id=50,41 # bp-sh ma repl and bp-sh discovery
## 0.8536 se=0.157
run_ldscore --rg --gwas-id=50,42 # bp-sh ma repl and bp-sh eur replication
## 1.3883 se=0.1
run_ldscore --rg --gwas-id=50,49 # bp-sh ma repl and bp-sh ma eur
## 1.0981 se=0.077


