# Aim of this script is to estimate heritability of Back pain after
# shared heredity subtraction and its genetic correlation with shared heredity

library(data.table)

setwd("/mnt/polyomica/projects/bp-sh/src/elgaeva_src/04_calculate_bp-sh/")

source("../../shared_heredity/00_core_functions/heritability_of_linear_combination.R")
source("../../shared_heredity/00_core_functions/gcor_a1_a2.R")
source("../../shared_heredity/00_core_functions/gcov_for_linear_comb_with_i_trait.R")

aa <- read.table('../../../data/01_sh//alphas.txt', row.names = 1)

alphas <- as.numeric(aa[2, ])

n_traits <- c(1:length(alphas))

phem <- read.table('../../../data/01_sh/pheno_corr_matrix.txt', row.names = 1,  check.names = F)

gcov <- read.table('../../../data/01_sh/gene_cov_matrix.txt', row.names = 1,  check.names = F)

slope <- sapply(n_traits, function(x) cov_gi_alpha(a = alphas, i = x, covm = as.matrix(gcov))/H2(alphas, covm = gcov, phem = phem))

position <- diag(length(alphas))

# Estimate heritability
tmp <- lapply(n_traits, function(x) H2(position[x, ] - alphas*slope[x], covm = gcov, phem = phem))
#[[1]]
#[1] 0.007463425 for hip pain-sh
#[[2]]
#[1] 0.01461107 for back pain-sh
#[[3]]
#[1] 0.006457976 for neck pain-sh
#[[4]]
#[1] 0.02428392 for knee pain-sh
#[[5]]
#[1] 0.02690786 for head pain-sh
#[[6]]
#[1] 0.005746272 for stomach pain-sh


# Estimate pairwise genetic correlations for traits-SH and SH
tmp2 <- lapply(n_traits, function(x) cor_gi_a1_a2(a1 = alphas, a2 = position[x, ] - alphas*slope[x], covm = gcov))
# [[1]]
# [1] -2.828575e-15 rg for sh and hip pain-sh
# [[2]]
# [1] -2.243392e-15 rg for sh and back pain-sh
# [[3]]
# [1] -7.030881e-16 rg for sh and neck pain-sh
# [[4]]
# [1] -3.766432e-15 rg for sh and knee pain-sh
#[[5]]
#[1] 1.99153e-15 rg for sh and head pain-sh
#[[6]]
#[1] 2.932129e-15 rg for sh and stomach pain-sh


