# Aim of this script is to estimate heritability of adjusted BP-SH
# and its genetic correlation with shared heredity, back pain, BP-SH, CKP and CHP 

library(data.table)

setwd("/mnt/polyomica/projects/bp-sh/src/elgaeva_src/13_calculate_bp-sh_adj/")

source("../../shared_heredity/00_core_functions/heritability_of_linear_combination.R")
source("../../shared_heredity/00_core_functions/gcor_a1_a2.R")
source("../../shared_heredity/00_core_functions/gcov_for_linear_comb_with_i_trait.R")

# coefficients for SH
alphas <- read.table('../../../data/01_sh/alphas.txt', row.names = 1)
alphas <- as.numeric(alphas[2, ])
# alphas <- c(alphas, c(0, 0))

# coefficients for adjusted BP-SH
aa <- c(1.0297036, 0.1735681, 0.1063936) # "UGC of CBP", "Knee", "Head" 
aa <- c(0, 0, 0, aa[2], aa[3], 0, 0, aa[1])

n_traits <- c(1:length(alphas))

phem <- read.table('../../../data/09_heatmaps/pheno_corr_matrix.txt', row.names = 1,  check.names = F) # "Hip", "Back", "Neck", "Knee", "Head", "Stomach", "SGC", "UGC of CBP"
gcov <- read.table('../../../data/09_heatmaps/gene_cov_matrix.txt', row.names = 1,  check.names = F) # "Hip", "Back", "Neck", "Knee", "Head", "Stomach", "SGC", "UGC of CBP"

position <- diag(length(alphas))

# Estimate heri tability
H2(aa, covm = gcov, phem = phem)
# [1] 0.01297899

# Estimate pairwise genetic correlations for adjusted BP-SH and other traits
tmp <- lapply(n_traits, function(x) cor_gi_a1_a2(a1 = position[x, ], a2 = aa, covm = gcov))
#[[1]]
#[1] 1.810107e-01 rg for adjusted bp-sh and hip pain
#[[2]]
#[1] 7.103524e-01 rg for adjusted bp-sh and back pain-sh
#[[3]]
#[1] 3.483217e-01 rg for adjusted bp-sh and neck pain-sh
#[[4]]
#[1] -4.873228e-08 rg for adjusted bp-sh and knee pain-sh
#[[5]]
#[1] 4.968217e-08 rg for adjusted bp-sh and head pain-sh
#[[6]]
#[1] 8.758968e-02 rg for adjusted bp-sh and stomach pain-sh
#[[7]]
#[1] 3.350627e-01 rg for adjusted bp-sh and sh
#[[7]]
#[1] 9.301499e-01 rg for adjusted bp-sh and bp-sh


