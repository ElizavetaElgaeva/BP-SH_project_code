# Aim of this script is to estimate genetic correlation between the
# shared heredity (meta-analysis of Europeans) and
# GIP1 (meta-analysis of Europeans) from our previous study (https://doi.org/10.1038/s42003-020-1051-9)
# using LD Score implementd in GWAS-Map (from test container)

# Estimate pairwise genetic correlation
run_ldscore --rg --gwas-id=46,47 
## 0.984 se=0.001

