# Aim of this script is to reformat input GWAS data on BP-SH Eur MA for COJO use

setwd("/mnt/polyomica/projects/bp-sh/data/")
bpsh <- data.table::fread("./00_upload_to_db/05_bp-sh_ma_eur/unification_results/BP-SH_ma_eur_output_done.csv")

bpsh <- bpsh[,c('rs_id','ea','ra','eaf','beta','se','p','n')]
colnames(bpsh) <- c('SNP','A1','A2','freq','b','se','p','N')

data.table::fwrite(bpsh, file = "./07_cojo/02_bp-sh/BP-SH_ma_eur_output_done.for_cojo", sep = '\t')

