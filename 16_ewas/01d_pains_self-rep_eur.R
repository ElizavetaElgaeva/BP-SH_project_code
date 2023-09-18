# Aim of this script is to prepare CBP phenotype 
# Self-reported Europeans only

library(data.table)
library(dplyr)

setwd("/home/lima/Zor/BACK_PAIN/New_pain_traits/Pheno/")

# Read CBP Phenotypic data
cbp <- fread("/home/lima/Zor/BACK_PAIN/New_pain_traits/Pheno/CBPw.phe", data.table = F)

# Read CKP Phenotypic data
ckp <- fread("/home/lima/Zor/BACK_PAIN/New_pain_traits/Pheno/CKPw.phe", data.table = F)

# Read CNP Phenotypic data
cnp <- fread("/home/lima/Zor/BACK_PAIN/New_pain_traits/Pheno/CNPw.phe", data.table = F)

# Read CHP Phenotypic data
chp <- fread("/home/lima/Zor/BACK_PAIN/New_pain_traits/Pheno/CHPw.phe", data.table = F)

table(is.na(cbp$chron))
table(is.na(ckp$chron))
table(is.na(cnp$chron))
table(is.na(chp$chron))

# Save data
cbp <- cbp[, c("ID", "ID", "chron")]
fwrite(cbp, file = "/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CBP_self-rep_eur_gcta_v2.txt", sep = " ",  col.names = F)

ckp <- ckp[, c("ID", "ID", "chron")]
fwrite(ckp, file = "/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CKP_self-rep_eur_gcta.txt", sep = " ",  col.names = F)

cnp <- cnp[, c("ID", "ID", "chron")]
fwrite(cnp, file = "/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CNP_self-rep_eur_gcta.txt", sep = " ",  col.names = F)

chp <- chp[, c("ID", "ID", "chron")]
fwrite(chp, file = "/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CHP_self-rep_eur_gcta.txt", sep = " ",  col.names = F)

# Check prevalence
covar <- fread("/home/common/projects/varicose_project/2023/pr59345_bcovar_sex_batch.txt")
qcovar <- fread("/home/common/projects/varicose_project/2023/pr59345_qcovar_pc_age.txt")
exomes <- fread("/home/common/DataStorage/UKBB/Project_59345/Exome/ukb23155/white_poly_qc/ukb23155_w_p_qc.psam")
imp <- fread("/home/common/DataStorage/UKBB/Project_59345/Imputed/PGEN/ukb_imp_chr10_v3.psam")
sample <- fread("/home/common/DataStorage/UKBB/Project_59345/Imputed/BGEN-1.3_hg38/ukb_imp_v3_hg38.sample")
sparse <- fread("/home/common/DataStorage/UKBB/Project_59345/GRM/ukb_ea_sp.grm_pr59345.ids")

# CBP
# add common filters
i <- which(cbp$ID %in% covar$V1)
cbp <- cbp[i, ]
table(cbp$ID %in% qcovar$V1) # TRUE
sp <- which(cbp$ID %in% sparse$V1)
cbp <- cbp[sp, ]

# exomes
e <- which(cbp$ID %in% exomes$IID)
table(cbp[e, c("chron")])
#     0      1
#148708  31415
31415/(31415+148708) # 0.1744086
# 180123 N total

# imputed
s <- which(cbp$ID %in% sample$"ID_1")
cbp <- cbp[s, ]
g <- which(cbp$ID %in% imp$IID)
table(cbp[g, c("chron")])
#     0      1
#362318  77953
77953/(77953+362318) # 0.1770569
# 440271 total


# CKP
# add common filters
i <- which(ckp$ID %in% covar$V1)
ckp <- ckp[i, ]
table(ckp$ID %in% qcovar$V1) # TRUE
sp <- which(ckp$ID %in% sparse$V1)
ckp <- ckp[sp, ]

# exomes
e <- which(ckp$ID %in% exomes$IID)
table(ckp[e, c("chron")])
#     0      1
#149984  30221
30221/(30221+149984) # 0.1677034
# 180205 N total

# imputed
s <- which(ckp$ID %in% sample$"ID_1")
ckp <- ckp[s, ]
g <- which(ckp$ID %in% imp$IID)
table(ckp[g, c("chron")])
#     0      1
#365662  74826
74826/(74826+365662) # 0.1698707
# 440488 N total

# CNP
# add common filters
i <- which(cnp$ID %in% covar$V1)
cnp <- cnp[i, ]
table(cnp$ID %in% qcovar$V1) # TRUE
sp <- which(cnp$ID %in% sparse$V1)
cnp <- cnp[sp, ]

# exomes
e <- which(cnp$ID %in% exomes$IID)
table(cnp[e, c("chron")])
#     0      1
#151726  28383
28383/(28383+151726) # 0.1575879
# 180109 N total

# imputed
s <- which(cnp$ID %in% sample$"ID_1")
cnp <- cnp[s, ]
g <- which(cnp$ID %in% imp$IID)
table(cnp[g, c("chron")])
#     0      1
#369609  70598
70598/(70598+369609) # 0.1603746
# 440207 N total

# CHP
# add common filters
i <- which(chp$ID %in% covar$V1)
chp <- chp[i, ]
table(chp$ID %in% qcovar$V1) # TRUE
sp <- which(chp$ID %in% sparse$V1)
chp <- chp[sp, ]

# exomes
e <- which(chp$ID %in% exomes$IID)
table(chp[e, c("chron")])
#     0      1
#164590  15717
15717/(15717+164590) # 0.087168
# 180307 N total

# imputed
s <- which(chp$ID %in% sample$"ID_1")
chp <- chp[s, ]
g <- which(chp$ID %in% imp$IID)
table(chp[g, c("chron")])
#     0      1
#401515  39207
39207/(39207+401515) # 0.08896084
# 440722 N total


# Calculate phenotypic correlations

nonrel <- fread("/home/common/DataStorage/UKBB/Project_59345/basket_2006911/2022/21000_22021/pr59345_white_unrel.lst")
cbp <- fread("/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CBP_self-rep_eur_gcta_v2.txt",  data.table = F)





colnames(cbp) <- c("ID", "ID", "chron")
colnames(ckp) <- c("ID", "ID", "chron")
colnames(chp) <- c("ID", "ID", "chron")
colnames(cnp) <- c("ID", "ID", "chron")

# Add commmon filters 
# ....
n <- which(cbp$ID %in% nonrel$V1)
cbp <- cbp[n, ]

# Exomes

cbp_ckp_e <- inner_join(cbp_e[, 2:3], ckp_e[, 2:3])
cor(cbp_ckp_e$CBP, cbp_ckp_e$CKP)
0.1995229

cbp_cnp_e <- inner_join(cbp_e[, 2:3], cnp_e[, 2:3])
cor(cbp_cnp_e$CBP, cbp_cnp_e$CNP)
0.2652037

cbp_chp_e <- inner_join(cbp_e[, 2:3], chp_e[, 2:3])
cor(cbp_chp_e$CBP, cbp_chp_e$CHP)
0.2488513

ckp_cnp_e <- inner_join(ckp_e[, 2:3], cnp_e[, 2:3])
cor(ckp_cnp_e$CKP, ckp_cnp_e$CNP)
0.1752323

ckp_chp_e <- inner_join(ckp_e[, 2:3], chp_e[, 2:3])
cor(ckp_chp_e$CKP, ckp_chp_e$CHP)
0.2318095

cnp_chp_e <- inner_join(cnp_e[, 2:3], chp_e[, 2:3])
cor(cnp_chp_e$CNP, cnp_chp_e$CHP)
0.1815975

