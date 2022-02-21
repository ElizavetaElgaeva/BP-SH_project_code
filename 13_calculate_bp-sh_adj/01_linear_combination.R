# This script is to adjust UGIT CBP for chronic headache and knee pain in Discovery study

library(data.table)
library(dplyr)

setwd('/mnt/polyomica/projects/bp-sh/src/elgaeva_src/13_calculate_bp-sh_adj')

source("../../shared_heredity/00_core_functions/linear_combination_v3.R")
source("../../shared_heredity/00_core_functions/gcov_for_linear_comb_with_i_trait.R")
source("../../shared_heredity/00_core_functions/heritability_of_linear_combination.R")

gwas_files <- c('../../../data/00_upload_to_db/00_discovery/02_bp-sh/unification_results/BP-SH_disc_output_done.csv',
		'../../../data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Knee/Knee_output_done.csv',
		'../../../data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Head/Head_output_done.csv')

gwas <- lapply(gwas_files, fread)

aa <- c(1.0297036, 0.1735681, 0.1063936) # "UGC of CBP", "Knee", "Head"
phem <- read.table('../../../data/09_heatmaps/pheno_corr_matrix.txt', row.names = 1,  check.names = F) # "Hip", "Back", "Neck", "Knee", "Head", "Stomach", "SGC", "UGC of CBP"
colnames(phem) <- rownames(phem) <- c("Hip", "Back", "Neck", "Knee", "Head", "Stomach", "SGC", "UGC of CBP")
phem <- phem[c("UGC of CBP", "Knee", "Head"), c("UGC of CBP", "Knee", "Head")]
gcov <- read.table('../../../data/09_heatmaps/gene_cov_matrix.txt', row.names = 1,  check.names = F) # "Hip", "Back", "Neck", "Knee", "Head", "Stomach", "SGC", "UGC of CBP"
colnames(gcov) <- rownames(gcov) <- c("Hip", "Back", "Neck", "Knee", "Head", "Stomach", "SGC", "UGC of CBP")
gcov <- gcov[c("UGC of CBP", "Knee", "Head"), c("UGC of CBP", "Knee", "Head")]

rs_id <- lapply(gwas, function(x) x$rs_id)
snps <- Reduce(intersect, rs_id)
ind <- lapply(rs_id, function(x) match(snps, x))

gwas_reordered <- lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]], ])
z <- sapply(gwas_reordered, function(x) x$z)
eaf <- gwas_reordered[[1]]$eaf
sample_size <- sapply(gwas_reordered, function(x) x$n)

n_traits <- c(1:length(aa))

tr_sh_adj_gwas <- lapply(n_traits, function(x) GWAS_linear_combination_Z_based(a = aa, Z = z, covm = as.matrix(phem), N = sample_size, eaf = eaf))
tr_sh_adj_gwas <- lapply(tr_sh_adj_gwas, function(x) mutate(x, Z = b/se, p = pchisq(Z^2, 1, low = F)))
tr_sh_adj_gwas <- lapply(tr_sh_adj_gwas, function(x) mutate(x, SNP = gwas_reordered[[1]]$rs_id))
tr_sh_adj_gwas <- lapply(tr_sh_adj_gwas, function(x) mutate(x, A1 = gwas_reordered[[1]]$ea, A2 = gwas_reordered[[1]]$ra, chr = gwas_reordered[[1]]$chr, pos = gwas_reordered[[1]]$bp, eaf = gwas_reordered[[1]]$eaf))

head(tr_sh_adj_gwas[[2]], n = 2)


fwrite(tr_sh_adj_gwas[[2]], 
       row.names = F,
       file = '/mnt/polyomica/projects/bp-sh/data/13_bp-sh_adj/bp-sh_adj_discovery/BP-SH_adj_discovery.txt',
       sep = '\t')

