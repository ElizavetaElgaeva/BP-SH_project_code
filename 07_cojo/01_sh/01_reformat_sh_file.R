# Aim of this script is to reformat input GWAS data on SH for COJO use

setwd("/mnt/polyomica/projects/bp-sh/data/")

# Eur MA
sh <- data.table::fread("./00_upload_to_db/03_sh_ma_eur/unification_results/SH_ma_eur_output_done.csv")

sh <- sh[, c('rs_id','ea','ra','eaf','beta','se','p','n')]
colnames(sh) <- c('SNP','A1','A2','freq','b','se','p','N')

data.table::fwrite(sh, file = "./07_cojo/01_sh/SH_ma_eur_output_done.for_cojo", sep = '\t')

# Discovery
sh <- data.table::fread("./00_upload_to_db/00_discovery/01_sh/unification_results/SH_disc_output_done.csv")

sh <- sh[, c('rs_id','ea','ra','eaf','beta','se','p','n')]
colnames(sh) <- c('SNP','A1','A2','freq','b','se','p','N')

data.table::fwrite(sh, file = "./07_cojo/01_sh/SH_disc_output_done.for_cojo", sep = '\t')

