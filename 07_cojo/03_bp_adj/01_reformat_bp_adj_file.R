# Aim of this script is to reformat input GWAS data on BP adjusted for COJO use

setwd("/mnt/polyomica/projects/bp-sh/data/")

# Discovery
bpadj <- data.table::fread("./00_upload_to_db/00_discovery/03_bp_adj/unification_results/BP_adj_disc_output_done.csv")

bpadj <- bpadj[, c('rs_id','ea','ra','eaf','beta','se','p','n')]
colnames(bpadj) <- c('SNP','A1','A2','freq','b','se','p','N')

data.table::fwrite(bpadj, file = "./07_cojo/03_bp_adj/BP_adj_disc_output_done.for_cojo", sep = '\t')

