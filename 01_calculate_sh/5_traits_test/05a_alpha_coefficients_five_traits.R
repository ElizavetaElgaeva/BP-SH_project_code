# Aim of this script is to estimate alpha coefficients in linear combination
# for five traits and estimate the expected genetic correlation between Shared Heredity
# obtained for five and six pain traits

source("../../shared_heredity/00_core_functions/shared_heredity.R")
source("../../shared_heredity/00_core_functions/gcor_a1_a2.R")

path <- '/mnt/polyomica/projects/bp-sh/data/01_sh/' # path to the data

CorPhenTr <- as.matrix(read.table(paste0(path, 'pheno_corr_matrix.txt'), check.names = F)) # load matrix of phenotypic correlations
CorPhenTr <- CorPhenTr[-5, -5] # exclude Head pain
A0 <- as.matrix(read.table(paste0(path, 'gene_cov_matrix.txt'), check.names = F)) # load matrix of genetic covariance
A1 <- A0[-5, -5] # exclude Head pain

# Estimate alpha coefficients for five traits
x <- shared_heredity(CovGenTr = A1, CorPhenTr = CorPhenTr)
x$alphas
alphas5 <- as.numeric(x$alphas[2, ])
alphas5
# 0.2837482 0.4341831 0.3981446 0.3412508 0.2000272
alphas5 <- append(alphas5, 0.0000000, 4) # add zero coeffitient for Head pain

# Estimate the expected genetic correlation between Shared Heredity obtained for five and six pain traits
aa <- read.table(paste0(path, 'alphas.txt'), row.names = 1) # load the linear combination coeffitients for six pain traits
alphas6 <- as.numeric(aa[2, ])

cor_gi_a1_a2(a1 = alphas5, a2 = alphas6, covm = A0)
# 0.9915223
