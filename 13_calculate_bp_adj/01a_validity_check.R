# Aim of this script is to estimate heritability of adjusted BP-SH
# and its genetic correlation with shared heredity, back pain, BP-SH, CKP and CHP 

library(data.table)

setwd("/mnt/polyomica/projects/bp-sh/src/elgaeva_src/13_calculate_bp_adj/")

source("../../shared_heredity/00_core_functions/heritability_of_linear_combination.R")
source("../../shared_heredity/00_core_functions/gcor_a1_a2.R")
source("../../shared_heredity/00_core_functions/gcov_for_linear_comb_with_i_trait.R")
source("../../shared_heredity/00_core_functions/cor_g_a.R")

# coefficients for SH
alphas <- read.table('../../../data/01_sh/alphas.txt', row.names = 1)
alphas <- as.numeric(alphas[2, ])

# coefficients for adjusted BP-SH
aa <- c(-0.071018123, 0.854891273, -0.790903139, 0.017435276, 0.001350986, -0.067428072) # "Hip", "Back", "Neck", "Knee", "Head", "Stomach"

n_traits <- c(1:length(alphas))

phem <- read.table('../../../data/01_sh/pheno_corr_matrix.txt', row.names = 1,  check.names = F) # "Hip", "Back", "Neck", "Knee", "Head", "Stomach"
gcov <- read.table('../../../data/01_sh/gene_cov_matrix.txt', row.names = 1,  check.names = F) # "Hip", "Back", "Neck", "Knee", "Head", "Stomach"

slope <- sapply(n_traits, function(x) cov_gi_alpha(a = alphas, i = x, covm = as.matrix(gcov))/H2(alphas, covm = gcov, phem = phem))

position <- diag(length(alphas))

# Estimate heritability
H2(aa, covm = gcov, phem = phem)
# [1] 0.007168787

# Estimate pairwise genetic correlations for adjusted BP-SH and other traits
cor_gi_a1_a2(a1 = alphas, a2 = aa, covm = gcov[c(1:6), c(1:6)]) 
#[1] 0.151285 rg for adjusted bp-sh and sh
cor_gi_a1_a2(a1 = position[2, ] - alphas*slope[2], a2 = aa, covm = gcov[c(1:6), c(1:6)])
#[1] 0.8021326 rg for adjusted bp-sh and bp-sh
tmp <- lapply(n_traits, function(x) cor_gi_alfa(i = n_traits[x], a = aa, covm = gcov))
#[[1]]
#[1] 3.604807e-11 rg for adjusted bp-sh and hip pain
#[[2]]
#[1] 0.4823323 rg for adjusted bp-sh and back pain
#[[3]]
#[1] -2.194532e-10 rg for adjusted bp-sh and neck pain
#[[4]]
#[1] 1.845365e-10 rg for adjusted bp-sh and knee pain
#[[5]]
#[1] 2.357728e-10 rg for adjusted bp-sh and head pain
#[[6]]
#[1] 7.117874e-11 rg for adjusted bp-sh and stomach pain

