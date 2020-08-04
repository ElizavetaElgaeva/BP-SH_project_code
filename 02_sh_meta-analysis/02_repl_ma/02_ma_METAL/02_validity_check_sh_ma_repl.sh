# Aim of this script is to estimate heritability of
# Replication meta-analysis of shared heredity and its genetic correlations with
# Shared heredity from Discovery study, European Replication study and European meta-analysis
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for SH MA Repl
run_ldscore --h2 --gwas-id=48
## h2=0.0731  se=0.0037

# Estimate pairwise genetic correlations for SH MA Repl
run_ldscore --rg --gwas-id=48,19 # sh ma repl and sh discovery
## 0.9783 se=0.03
run_ldscore --rg --gwas-id=48,26 # sh ma repl and sh eur replication
## 1.0757 se=0.005
run_ldscore --rg --gwas-id=48,47 # sh ma repl and sh ma eur
## 1.0228 se=0.012

