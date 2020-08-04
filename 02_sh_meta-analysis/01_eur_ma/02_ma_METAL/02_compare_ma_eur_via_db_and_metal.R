# Aim of this script is to compare the results of meta-analysis of Europeans
# obtained using GWAS-MAP and METAL

library(data.table)
library(dplyr)

setwd('/mnt/polyomica/projects/bp-sh/src/elgaeva_src/02_sh_meta-analysis/01_eur_ma/02_ma_METAL')

# Load the meta-analysis results
ma_db <- fread('/mnt/polyomica/projects/bp-sh/data/02_sh_ma/01_eur_ma/01_ma_via_db/SH_EUR_MA.csv', data.table = FALSE)
colnames(ma_db)
dim(ma_db)

ma_metal <- fread('/mnt/polyomica/projects/bp-sh/data/02_sh_ma/01_eur_ma/02_ma_METAL/SH_ma_eur_metal.txt', data.table = FALSE)
colnames(ma_metal)
dim(ma_metal)

# Reorder the data
table(ma_db$rs_id %in% ma_metal$MarkerName)
ind <- match(ma_db$rs_id, ma_metal$MarkerName)
ma_metal <- ma_metal[ind, ]

# Create a subset
ind <- sample(1:nrow(ma_db), 100000, replace = FALSE)

# Flip the alleles
for(i in ind){
	if((ma_db$ra[i] == toupper(ma_metal$Allele1[i])) & (ma_db$ea[i] == toupper(ma_metal$Allele2[i]))){
		ma_metal$Allele1[i] <- tolower(ma_db$ea[i])
		ma_metal$Allele2[i] <- tolower(ma_db$ra[i])
		ma_metal$Effect[i] <- ma_metal$Effect[i]*(-1)
		ma_metal$Freq1[i] <- 1 - ma_metal$Freq1[i]
	}
	if((ma_db$ea[i] != toupper(ma_metal$Allele1[i])) | (ma_db$ra[i] != toupper(ma_metal$Allele2[i]))){
		print(paste0('Check for triallelic SNP in ', ma_metal$MarkerName[i], '\n'))
	}
}

# Compare SNP effects
out <- '/mnt/polyomica/projects/bp-sh/data/02_sh_ma/01_eur_ma/02_ma_METAL/beta_plot.png'
png(out, height = 720, width= 720)
plot(ma_db$beta_ma[ind],
     ma_metal$Effect[ind],
     xlab = 'Beta from GWAS-MAP',
     ylab = 'Beta from METAL',
     main = 'SH MA EUR betas')
dev.off()

# Compare SNP SE of the effects
out <- '/mnt/polyomica/projects/bp-sh/data/02_sh_ma/01_eur_ma/02_ma_METAL/se_plot.png'
png(out, height = 720, width= 720)
plot(ma_db$se_ma[ind],
     ma_metal$StdErr[ind],
     xlab = 'SE from GWAS-MAP',
     ylab = 'SE from METAL',
     main = 'SH MA EUR standard errors')
dev.off()

# Compare SNP p-values
out <- '/mnt/polyomica/projects/bp-sh/data/02_sh_ma/01_eur_ma/02_ma_METAL/p_plot.png'
png(out, height = 720, width= 720)
plot(ma_db$p_ma[ind],
     ma_metal$"P-value"[ind],
     xlab = 'P-value from GWAS-MAP',
     ylab = 'P-value from METAL',
     main = 'SH MA EUR p-values')
dev.off()

