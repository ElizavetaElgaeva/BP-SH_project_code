# Aim of this script is to create a list of SNPs from SH GWAS to perform SMR-HEIDI analysis

library(data.table)

setwd('/mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi/')

t <- fread('../../10_tables/st1_sh_extended.tsv', data.table = F)
t <- t[, c('SNP', 'CHR', 'POS')]
not_repl <- match(c('rs34802737', 'rs10137555'), t$SNP)
t <- t[-not_repl, ]
colnames(t) <- c('rs_id', 'chr', 'bp')
fwrite(t, file = './snp_list.csv', sep = ',')
