# Aim of this script is to prepare CBP phenotype 
# Self-reported Europeans only

library(data.table)
library(dplyr)

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/")

# Read list of White individuals (self-reported)
new_white <- fread("/home/common/DataStorage/UKBB/Project_18219/addons/pr18219_white.lst", data.table = F)

# Read Phenotypic data
pheno <- fread("/home/freydin/Pheno/mv_mar18_chron.phen.txt", data.table = F)

# Match the data
table(new_white$V1 %in% pheno$IID)
# TRUE
# 459161

final <- inner_join(pheno, new_white, by = c("IID" = "V1"))
dim(final)

table(is.na(final$back))
# FALSE   TRUE
#447466  11695

final <- final[!is.na(final$back), ]
dim(final)
#[1] 447466     26

table(final$back)
#     0      1
#367058  80408

80408/447466
# 0.1796963 prevalence

# Save data
readme_final <- "UKBB table with CBP phenotype. Self-reported Europeans only, N=447K"
save(final, readme_final, file = "/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/chronic_self-rep_eur_pheno.RData")

cbp <- final[, c("IID", "IID", "back")]
fwrite(cbp, file = "/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CBP_self-rep_eur_gcta.txt", sep = " ",  col.names = F)

sex_batch <-  final[, c("IID", "IID", "Sex", "batch")]
fwrite(sex_batch, file = "/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CBP_bcovar_sex_batch_self-rep_eur_gcta.txt", sep = " ",  col.names = F)

pcs <- final[ c("IID", "IID", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", "Age")]
fwrite(pcs, file = "/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CBP_qcovar_pc1_10_age_self-rep_eur_gcta.txt", sep = " ",  col.names = F)

