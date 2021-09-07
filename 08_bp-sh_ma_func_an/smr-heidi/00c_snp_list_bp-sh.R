# Aim of this script is to create a list of SNPs from BP-SH GWAS to perform SMR-HEIDI analysis

library(data.table)

setwd('/mnt/polyomica/projects/bp-sh/data/08_bp-sh_ma_func_an/smr-heidi/')

t <- fread('../../10_tables/st1_bp-sh_extended.tsv', data.table = F)
t <- t[, c('SNP', 'CHR', 'POS')]
colnames(t) <- c('rs_id', 'chr', 'bp')
fwrite(t, file = './snp_list.csv', sep = ',')
