# Aim of this script is to estimate alpha coefficients in linear combination
# for each trait

source("/mnt/polyomica/projects/bp-sh/src/shared_heredity/00_core_functions/shared_heredity.R")

setwd("/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/")

load("four_pains_GIPs.RData")

CorPhenTr <- as.matrix(four_pains$cor_y) # load matrix of phenotypic correlations
A0 <- as.matrix(four_pains$cov_g) # load matrix of genetic covariance

# Reorder and exclude unnecessry column 
CorPhenTr <- CorPhenTr[c('back', 'neck', 'knee', 'hip'), c('back', 'neck', 'knee', 'hip')]
A0 <- A0[c('back', 'neck', 'knee', 'hip'), c('back', 'neck', 'knee', 'hip')]

# Set names for output files
output_alpha <- 'alphas_4_pains.txt'
output_w <- 'weights_4_pains.txt'
#output_gip <- 'GIP1.txt'

# Estimate alpha coefficients
x <- shared_heredity(CovGenTr = A0, CorPhenTr = CorPhenTr)
print(x)  
write.table(x$alphas, output_alpha, quote = F)
write.table(x$weights, output_w, quote = F)
#write.table(x$GIPs$GIP_coeff[, 'GIP1'], output_gip, quote = F)
x$alphas
