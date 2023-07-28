# Aim of this script is to estimate heritability of Back pain after
# shared heredity subtraction and its genetic correlation with shared heredity
# and back pain

library(data.table)

setwd("/mnt/polyomica/projects/bp-sh/src/elgaeva_src/04_calculate_bp-sh/5_traits_test/")

source("../../../shared_heredity/00_core_functions/heritability_of_linear_combination.R")
source("../../../shared_heredity/00_core_functions/gcor_a1_a2.R")
source("../../../shared_heredity/00_core_functions/gcov_for_linear_comb_with_i_trait.R")

aa <- read.table('../../../../data/01_sh/5_traits_test/alphas.txt', row.names = 1)

alphas <- as.numeric(aa[2, ])
alphas <- append(x=alphas, values=0, after=1)

n_traits <- c(1:length(alphas))

phem <- read.table('../../../../data/01_sh/pheno_corr_matrix.txt', row.names = 1,  check.names = F)

gcov <- read.table('../../../../data/01_sh/gene_cov_matrix.txt', row.names = 1,  check.names = F)

slope <- sapply(n_traits, function(x) cov_gi_alpha(a = alphas, i = x, covm = as.matrix(gcov))/H2(alphas, covm = gcov, phem = phem))

position <- diag(length(alphas))

# Estimate heritability
tmp <- lapply(n_traits, function(x) H2(position[x, ] - alphas*slope[x], covm = gcov, phem = phem))
#[[1]]
#[1] 0.007568529 for hip pain-sh
#[[2]]
#[1] 0.01531052 for back pain-sh
#[[3]]
#[1] 0.01010604 for neck pain-sh
#[[4]]
#[1] 0.02112433 for knee pain-sh
#[[5]]
#[1] 0.02675848 for head pain-sh
#[[6]]
#[1] 0.00548701 for stomach pain-sh


# Estimate pairwise genetic correlations for traits-SH and SH
tmp2 <- lapply(n_traits, function(x) cor_gi_a1_a2(a1 = alphas, a2 = position[x, ] - alphas*slope[x], covm = gcov))
# [[1]]
# [1] -6.239545e-15 rg for sh and hip pain-sh
# [[2]]
# [1] -5.760481e-15 rg for sh and back pain-sh
# [[3]]
# [1] -7.776451e-15 rg for sh and neck pain-sh
# [[4]]
# [1] -1.020399e-14 rg for sh and knee pain-sh
#[[5]]
#[1] -7.711336e-16 rg for sh and head pain-sh
#[[6]]
#[1] -3.360199e-15 rg for sh and stomach pain-sh

# Estimate pairwise genetic correlations for traits-SH and original pain traits
tmp3 <- lapply(n_traits, function(x) cor_gi_a1_a2(a1 = position[x, ], a2 = position[x, ] - alphas*slope[x], covm = gcov))
#[[1]]
#[1] 0.4583455 rg for hip pain-sh and hip pain
#[[2]]
#[1] 0.5948529 rg for back pain and back pain-sh
#[[3]]
#[1] 0.4214898 rg for neck pain and neck pain-sh
#[[4]]
#[1] 0.5493585 rg for knee pain and knee pain-sh
#[[5]]
#[1] 0.8246766 rg for head pain and head pain-sh
#[[6]]
#[1] 0.547719 rg for stomach pain and stomach pain-sh


