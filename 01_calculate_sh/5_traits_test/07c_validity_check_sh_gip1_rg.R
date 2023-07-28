# Aim of this script is to estimate the expected genetic correlation between Shared Heredity
# obtained for five pain traits and GIP1 obtained for four musculosceletal pain traits

source("../../../shared_heredity/00_core_functions/gcor_a1_a2.R")

path <- '/mnt/polyomica/projects/bp-sh/data/01_sh/5_traits_test/' # path to the data

A0 <- as.matrix(read.table(paste0(path, 'gene_cov_matrix.txt'), check.names = F)) # load matrix of genetic covariance
aa <- read.table(paste0(path, 'alphas.txt'), row.names = 1) # load the linear combination coeffitients for five pain traits
alphas <- as.numeric(aa[2, ])

load('/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/20181010_GPCs.RData') # load data for GIP of four musculosceletal pain traits
ls()
four_pains$eigens # GIP coeffitients for four traits
gip1 <- as.numeric(four_pains$eigens[,1])
gip1 <- append(gip1, gip1[4], 0) # reordering the elements
gip1[5] <- 0.0000000 # make the length of the vector as the alphas vector

# Estimate the expected genetic correlation between Shared Heredity obtained for five pain traits and GIP1 of four pain traits
cor_gi_a1_a2(a1 = alphas, a2 = gip1, covm = A0)
# 0.9939956
