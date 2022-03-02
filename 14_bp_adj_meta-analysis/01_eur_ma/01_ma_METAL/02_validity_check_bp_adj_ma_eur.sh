# Aim of this script is to estimate heritability of
# European meta-analysis of back pain adjusted for all other original pain traits
# and calculate its genetic correlations with
# back pain adjusted for all other original pain traits from Discovery study and European Replication study
# using LD Score implementd in GWAS-Map (from test container)

# Estimate heritability for BP adjusted MA EUR
run_ldscore --h2 --gwas-id=55
## h2=0.0082  se=0.0011

# Estimate pairwise genetic correlations for BP adjusted MA EUR
run_ldscore --rg --gwas-id=55,53 # BP adjusted ma eur and BP adjusted discovery
## 0.945 se=0.043
run_ldscore --rg --gwas-id=55,54 # BP adjusted ma eur and BP adjusted eur replication
## 0.9311 se=0.056

