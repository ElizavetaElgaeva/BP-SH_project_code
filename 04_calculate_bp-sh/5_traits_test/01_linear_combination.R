# This script is to substract shared genetic component from Back pain in Discovery study

library(data.table)
library(dplyr)

setwd('/mnt/polyomica/projects/bp-sh/src/elgaeva_src/04_calculate_bp-sh/5_traits_test')

source("../../../shared_heredity/00_core_functions/linear_combination_v3.R")
source("../../../shared_heredity/00_core_functions/gcov_for_linear_comb_with_i_trait.R")
source("../../../shared_heredity/00_core_functions/heritability_of_linear_combination.R")

gwas_files <- c('../../../../data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Hip/Hip_output_done.csv',
		'../../../../data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Back/Back_output_done.csv',
		'../../../../data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Neck/Neck_output_done.csv',
		'../../../../data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Knee/Knee_output_done.csv',
		'../../../../data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Head/Head_output_done.csv',
		'../../../../data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Stom/Stom_output_done.csv')

gwas <- lapply(gwas_files, fread)

aa <- read.table('../../../../data/01_sh/5_traits_test/alphas.txt', row.names = 1)
aa$Alpha14 <- c(0, 0)
aa <- aa[, c('Alpha.13', 'Alpha14', 'Alpha.15', 'Alpha.16', 'Alpha.17', 'Alpha.18')]

phem <- read.table('../../../../data/01_sh/pheno_corr_matrix.txt', row.names = 1,  check.names = F)

gcov <- read.table('../../../../data/01_sh/gene_cov_matrix.txt', row.names = 1,  check.names = F)

rs_id <- lapply(gwas, function(x) x$rs_id)
snps <- Reduce(intersect, rs_id)
ind <- lapply(rs_id, function(x) match(snps, x))

gwas_reordered <- lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]], ])
z <- sapply(gwas_reordered, function(x) x$z)
eaf <- gwas_reordered[[1]]$eaf
sample_size <- sapply(gwas_reordered, function(x) x$n)


alphas <- as.numeric(aa[2, ])
n_traits <- c(1:length(alphas))
slope <- sapply(n_traits, function(x) cov_gi_alpha(a = alphas, i = x, covm = as.matrix(gcov))/H2(alphas, covm = gcov, phem = phem))
position <- diag(length(alphas))


tr_sh_gwas <- lapply(n_traits, function(x) GWAS_linear_combination_Z_based(a = position[x, ] - alphas*slope[x], Z = z, covm = as.matrix(phem), N = sample_size, eaf = eaf))
tr_sh_gwas <- lapply(tr_sh_gwas, function(x) mutate(x, Z = b/se, p = pchisq(Z^2, 1, low = F)))
tr_sh_gwas <- lapply(tr_sh_gwas, function(x) mutate(x, SNP = gwas_reordered[[1]]$rs_id))
tr_sh_gwas <- lapply(tr_sh_gwas, function(x) mutate(x, A1 = gwas_reordered[[1]]$ea, A2 = gwas_reordered[[1]]$ra, chr = gwas_reordered[[1]]$chr, pos = gwas_reordered[[1]]$bp, eaf = gwas_reordered[[1]]$eaf))

head(tr_sh_gwas[[2]], n = 2)


fwrite(tr_sh_gwas[[2]], 
       row.names = F,
       file = '/mnt/polyomica/projects/bp-sh/data/04_bp-sh/bp-sh_discovery/5_traits_test/BP-SH_discovery.txt',
       sep = '\t')

