# Aim of this script is to make final lists of eQTL GWASes for SMR-HEIDI analysis by locus

library(data.table)

setwd('/mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi')

gwas_files <- list.files(getwd(), pattern="locus_chr")
gwas <- lapply(gwas_files, fread)
gwas_filtered <- lapply(gwas, unique)

# Rewrite using filtered files
lapply(1:length(gwas_filtered), function(x) fwrite(gwas_filtered[[x]], row.names = F, col.names = F, file = gwas_files[[x]], sep = '\t'))
