# Aim of this script is to make a list of GWAS ids to build a full squared matrix of rg
# SH and BP-SH traits will be included
# only those traits, that have statistically significant (p-value < 3.424658e-05) and large enough (|rg| > 0.25) rg with SH or BP-SH will be selected

setwd('/mnt/polyomica/projects/bp-sh/data/')

library(data.table)

thr <- 0.05/(730*2) # significance threshold

rg <- fread("./03_sh_ma_func_an/gene_corrs/rg_sh_ma_eur.csv", data.table=F)
bprg <- fread("./08_bp-sh_ma_func_an/gene_corrs/rg_bp-sh_ma_eur.csv", data.table = F)

ids <- rg[rg$pval < 3.424658e-05 & abs(rg$rg) > 0.25, c("gwas_id_2")]
ids <- c(ids, bprg[bprg$pval < 3.424658e-05 & abs(bprg$rg) > 0.25, c("gwas_id_2")])
ids <- unique(ids)
length(ids)

ids <- c(ids, c(2269986, 2269987))
fwrite(as.data.table(ids), file = './03_sh_ma_func_an/gene_corrs/gwas_ids_for_sq_matrix.txt', row.names = F, col.names = F, sep = '\t')
