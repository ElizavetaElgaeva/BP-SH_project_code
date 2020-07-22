# Aim of this script is to estimate heritability of
# shared heredity and its genetic correlations with 6 pain traits

library(data.table)

setwd('/mnt/polyomica/projects/bp-sh/src/elgaeva_src/01_calculate_sh/')

source("../../shared_heredity/00_core_functions/heritability_of_linear_combination.R")
source("../../shared_heredity/00_core_functions/cor_g_a.R")

aa <- read.table('../../../data/01_sh/alphas.txt', row.names = 1)
phem <- read.table('../../../data/01_sh/pheno_corr_matrix.txt', row.names = 1,  check.names = F)
gcov <- read.table('../../../data/01_sh/gene_cov_matrix.txt', row.names = 1,  check.names = F)

alphas <- as.numeric(aa[2, ])

# Estimate heritability
H2(alphas, covm = gcov, phem = phem)
## 0.07024648

# Estimate pairwise genetic correlations for SH and 6 pain traits
cor_gi_alfa(a = alphas, i = 1, covm = gcov) # sh and Hip pain
## 0.8808927
cor_gi_alfa(a = alphas, i = 2, covm = gcov) # sh and Back pain
## 0.9022859
cor_gi_alfa(a = alphas, i = 3, covm = gcov) # sh and Neck pain
## 0.9387536
cor_gi_alfa(a = alphas, i = 4, covm = gcov) # sh and Knee pain
## 0.7824733
cor_gi_alfa(a = alphas, i = 5, covm = gcov) # sh and Head pain
## 0.5362005
cor_gi_alfa(a = alphas, i = 6, covm = gcov) # sh and Stomach pain
## 0.8197805



