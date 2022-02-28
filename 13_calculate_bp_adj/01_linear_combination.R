# This script is to adjust CBP for all other original chronic pain traits in Discovery study

library(data.table)
library(dplyr)

setwd('/mnt/polyomica/projects/bp-sh/src/elgaeva_src/13_calculate_bp_adj')

source("../../shared_heredity/00_core_functions/linear_combination_v3.R")
source("../../shared_heredity/00_core_functions/gcov_for_linear_comb_with_i_trait.R")
source("../../shared_heredity/00_core_functions/heritability_of_linear_combination.R")

gwas_files <- c('/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Hip/Hip_output_done.csv',
		'/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Back/Back_output_done.csv',
		'/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Neck/Neck_output_done.csv',
		'/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Knee/Knee_output_done.csv',
		'/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Head/Head_output_done.csv',
		'/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/unification_results/Stom/Stom_output_done.csv')

gwas <- lapply(gwas_files, fread)

aa <- c(0.854891273, -0.071018123, 0.017435276, 0.001350986, -0.067428072, -0.790903139) # "Back", "Hip", "Knee", "Head", "Stomach", "Neck" 
phem <- read.table('../../../data/01_sh/pheno_corr_matrix.txt', row.names = 1,  check.names = F) # "Hip", "Back", "Neck", "Knee", "Head", "Stomach"
colnames(phem) <- rownames(phem) <- c("Hip", "Back", "Neck", "Knee", "Head", "Stomach")
i <- match(colnames(phem), c("Back", "Hip", "Knee", "Head", "Stomach", "Neck"))
table(colnames(phem) == c("Back", "Hip", "Knee", "Head", "Stomach", "Neck")[i])
aa <- aa[i]
gcov <- read.table('../../../data/01_sh/gene_cov_matrix.txt', row.names = 1,  check.names = F)
colnames(gcov) <- rownames(gcov) <- colnames(phem)

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
       file = '/mnt/polyomica/projects/bp-sh/data/13_bp_adj/bp_adj_discovery/BP_adj_discovery.txt',
       sep = '\t')

