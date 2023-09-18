# Aim of this script is to obtain GWAS summary statistics for shared heredity in Discovery cohort

library(data.table)
library(dplyr)
source("/mnt/polyomica/projects/bp-sh/src/shared_heredity/00_core_functions/linear_combination_v3.R") 

setwd('/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/')
load("four_pains_GIPs.RData")
ls()

set1 <- c('/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CBP/CBP_self-rep_eur_exomes_set1.mac10_collaps_v2.fastGWA',
	  '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CNP/CNP_self-rep_eur_exomes_set1.mac10_collaps.fastGWA',
	  '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CKP/CKP_self-rep_eur_exomes_set1.mac10_collaps.fastGWA',
	  '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CHP/CHP_self-rep_eur_exomes_set1.mac10_collaps.fastGWA')

set2 <- c('/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CBP/CBP_self-rep_eur_exomes_set2.mac10_collaps_v2.fastGWA',
          '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CNP/CNP_self-rep_eur_exomes_set2.mac10_collaps.fastGWA',
          '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CKP/CKP_self-rep_eur_exomes_set2.mac10_collaps.fastGWA',
          '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CHP/CHP_self-rep_eur_exomes_set2.mac10_collaps.fastGWA')

set3 <- c('/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CBP/CBP_self-rep_eur_exomes_set3.mac10_collaps_v2.fastGWA',
          '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CNP/CNP_self-rep_eur_exomes_set3.mac10_collaps.fastGWA',
          '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CKP/CKP_self-rep_eur_exomes_set3.mac10_collaps.fastGWA',
          '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CHP/CHP_self-rep_eur_exomes_set3.mac10_collaps.fastGWA')

set4 <- c('/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CBP/CBP_self-rep_eur_exomes_set4.mac10_collaps_v2.fastGWA',
          '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CNP/CNP_self-rep_eur_exomes_set4.mac10_collaps.fastGWA',
          '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CKP/CKP_self-rep_eur_exomes_set4.mac10_collaps.fastGWA',
          '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CHP/CHP_self-rep_eur_exomes_set4.mac10_collaps.fastGWA')


raw_set <- c('/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CBP/CBP_self-rep_eur_exomes_v2.fastGWA',
	     '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CNP/CNP_self-rep_eur_exomes.fastGWA',
	     '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CKP/CKP_self-rep_eur_exomes.fastGWA',
	     '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CHP/CHP_self-rep_eur_exomes.fastGWA')


set1.maf0.01 <- c('/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CBP/CBP_self-rep_eur_exomes_set1.maf0.01_collaps_v2.fastGWA',
          '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CNP/CNP_self-rep_eur_exomes_set1.maf0.01_collaps.fastGWA',
          '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CKP/CKP_self-rep_eur_exomes_set1.maf0.01_collaps.fastGWA',
          '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CHP/CHP_self-rep_eur_exomes_set1.maf0.01_collaps.fastGWA')


gwas <- lapply(set1, fread)
gwas <- lapply(set2, fread)
gwas <- lapply(set3, fread)
gwas <- lapply(set4, fread)
gwas <- lapply(raw_set, fread)
gwas <- lapply(set1.maf0.01, fread)

aa <- read.table('alphas_4_pains.txt', row.names = 1)

covm <- four_pains$cor_y

# Reorder and exclude unnecessry column
covm <- covm[c('back', 'neck', 'knee', 'hip'), c('back', 'neck', 'knee', 'hip')]

rs_id <- lapply(gwas, function(x) x$"SNP")
snps <- Reduce(intersect, rs_id)
ind <- lapply(rs_id, function(x) match(snps, x))

gwas_reordered <- lapply(1:length(gwas), function(x) gwas[[x]][ind[[x]], ])
z <- sapply(gwas_reordered, function(x) x$"BETA"/x$"SE")
eaf <- gwas_reordered[[1]]$"AF1"
sample_size <- sapply(gwas_reordered, function(x) x$"N")

sh_gwas <- GWAS_linear_combination_Z_based(a = as.numeric(aa[2, ]), Z = z, covm = as.matrix(covm), N = sample_size, eaf = eaf)

sh_gwas <- mutate(sh_gwas, Z = b/se, p = pchisq(Z^2, 1, low = F))
sh_gwas <- mutate(sh_gwas, SNP = gwas_reordered[[1]]$"SNP")
sh_gwas <- mutate(sh_gwas, A1 = gwas_reordered[[1]]$"A1", A2 = gwas_reordered[[1]]$"A2", chr = gwas_reordered[[1]]$"CHR", pos = gwas_reordered[[1]]$"POS",
		  			 eaf = gwas_reordered[[1]]$"AF1")
# set1
head(sh_gwas, n = 2)
data.table::fwrite(
	sh_gwas, 
	file = '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/SGIT/SGIT_self-rep_eur_exomes_4_pains_set1.mac10_collaps.fastGWA',
	sep = '\t')

# set2
data.table::fwrite(
		   sh_gwas,
		   file = '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/SGIT/SGIT_self-rep_eur_exomes_4_pains_set2.mac10_collaps.fastGWA',
		   sep = '\t')

# set3
data.table::fwrite(
		   sh_gwas,
		   file = '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/SGIT/SGIT_self-rep_eur_exomes_4_pains_set3.mac10_collaps.fastGWA',
		   sep = '\t')

# set4
data.table::fwrite(
		   sh_gwas,
		   file = '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/SGIT/SGIT_self-rep_eur_exomes_4_pains_set4.mac10_collaps.fastGWA',
		   sep = '\t')

# raw_set
head(sh_gwas, n = 2)
data.table::fwrite(
		   sh_gwas,
		   file = '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/SGIT/SGIT_self-rep_eur_exomes_4_pains.fastGWA',
		   sep = '\t')

# set1.maf0.01
data.table::fwrite(
		   sh_gwas,
		   file = '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/SGIT/SGIT_self-rep_eur_exomes_4_pains_set1.maf0.01_collaps.fastGWA',
		   sep = '\t')


