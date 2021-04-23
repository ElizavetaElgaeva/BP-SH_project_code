# Aim of this script is to reformat input GWAS data on SH Eur MA for COJO use

setwd("/mnt/polyomica/projects/bp-sh/data/")
sh <- data.table::fread("./00_upload_to_db/03_sh_ma_eur/unification_results/SH_ma_eur_output_done.csv")

sh <- sh[,c('rs_id','ea','ra','eaf','beta','se','p','n')]
colnames(sh) <- c('SNP','A1','A2','freq','b','se','p','N')

data.table::fwrite(sh, file = "./07_cojo/01_sh/SH_ma_eur_output_done.for_cojo", sep = '\t')

