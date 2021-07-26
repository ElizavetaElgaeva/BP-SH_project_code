# Aim of this script is to create a full table representing all associated with SH loci
# and provide full information for inference that they are replicated

library(data.table)
library(dplyr)

setwd('/mnt/polyomica/projects/bp-sh/data')

# Creating an SNP-information-part of a table
backGWAS <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Back/Back_output_done.csv', data.table = F)
shGWAS <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/01_sh/unification_results/SH_disc_output_done.csv', data.table = F)
cojoFile <- fread(input = '/mnt/polyomica/projects/bp-sh/data/07_cojo/01_sh/sh_disc_all_chr_2.5e-8.jma.cojo', data.table = F)
st1 <- cojoFile$SNP
st1 <- as.data.frame(st1, stringsAsFactors=FALSE)
colnames(st1) <- c('SNP')
ind <- match(st1$SNP, shGWAS$rs_id)
shGWAS <- shGWAS[ind,]
st1 <- mutate(st1, CHR = cojoFile$Chr)
st1 <- mutate(st1, POS = cojoFile$bp)
st1 <- mutate(st1, REFERENCE_ALLELE = shGWAS$ra)
st1 <- mutate(st1, EFFECTIVE_ALLELE = shGWAS$ea)

# Creating a discovery-part of the table

ind <- match(st1$SNP, backGWAS$rs_id)
backGWAS <- backGWAS[ind,]
intercept <- 1.03758691779465 # LD Score intercept for SH Discovery calculated using --h2 option through the data base
st1 <- mutate(st1, MAF_DISC = pmin(shGWAS$eaf, (1 - shGWAS$eaf)))
st1 <- mutate(st1, EAF_DISC = shGWAS$eaf)
st1 <- mutate(st1, IMPUTATION_QUALITY_DISC = backGWAS$info)
st1 <- mutate(st1, BETA_DISC = shGWAS$beta)
st1 <- mutate(st1, SE_DISC = shGWAS$se)
st1 <- mutate(st1, PVAL_DISC = shGWAS$p)
st1 <- mutate(st1, PVAL_GC_DISC = pchisq(((st1$BETA_DISC / st1$SE_DISC)^2 / intercept), df = 1, lower.tail = F))
st1 <- mutate(st1, HWE_DISC = NA)
st1 <- mutate(st1, N_DISC = shGWAS$n)

# Creating a replication-EUR-part of the table

shEUR <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/EA/01_sh/unification_results/SH_ea_repl_output_done.csv', data.table = F)
ind <- match(st1$SNP, shEUR$rs_id)
shEUR <- shEUR[ind,]
backEUR <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/EA/00_original_traits/unification_results/Back/Back_output_done.csv', data.table = F)
ind <- match(st1$SNP, backEUR$rs_id)
backEUR <- backEUR[ind,]
st1 <- mutate(st1, BETA_EUR = shEUR$beta)
st1 <- mutate(st1, SE_EUR = shEUR$se)
st1 <- mutate(st1, PVAL_EUR = shEUR$p)
st1 <- mutate(st1, HWE_EUR = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_EUR = backEUR$info)
st1 <- mutate(st1, EAF_EUR = shEUR$eaf)
st1 <- mutate(st1, N_EUR = shEUR$n)

# Creating a replication-AFR-part of the table

shAFR <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/AA/01_sh/unification_results/SH_aa_repl_output_done.csv', data.table = F)
ind <- match(st1$SNP, shAFR$rs_id)
shAFR <- shAFR[ind,]
backAFR <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/AA/00_original_traits/unification_results/Back/Back_output_done.csv', data.table = F)
ind <- match(st1$SNP, backAFR$rs_id)
backAFR <- backAFR[ind,]
st1 <- mutate(st1, BETA_AFR = shAFR$beta)
st1 <- mutate(st1, SE_AFR = shAFR$se)
st1 <- mutate(st1, PVAL_AFR = shAFR$p)
st1 <- mutate(st1, HWE_AFR = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_AFR = backAFR$info)
st1 <- mutate(st1, EAF_AFR = shAFR$eaf)
st1 <- mutate(st1, N_AFR = shAFR$n)

# Creating a replication-ASI-part of the table

shASI <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/SA/01_sh/unification_results/SH_sa_repl_output_done.csv', data.table = F)
ind <- match(st1$SNP, shASI$rs_id)
shASI <- shASI[ind,]
backASI <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/SA/00_original_traits/unification_results/Back/Back_output_done.csv', data.table = F)
ind <- match(st1$SNP, backASI$rs_id)
backASI <- backASI[ind,]
st1 <- mutate(st1, BETA_ASI = shASI$beta)
st1 <- mutate(st1, SE_ASI = shASI$se)
st1 <- mutate(st1, PVAL_ASI = shASI$p)
st1 <- mutate(st1, HWE_ASI = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_ASI = backASI$info)
st1 <- mutate(st1, EAF_ASI = shASI$eaf)
st1 <- mutate(st1, N_ASI = shASI$n)

# Creating an MA-replication-part of the table

shMA <- fread(input = '/mnt/polyomica/projects/bp-sh/data/02_sh_ma/02_repl_ma/02_ma_METAL/SH_ma_repl_metal.txt')
tmp <- NULL
for (row in c(1:nrow(st1))) {
	if (nrow(shMA[MarkerName == st1$SNP[row]]) == 0) {
		tmp <- rbind(tmp, shMA[MarkerName == paste0(st1$CHR[row], ':', st1$POS[row], '_', st1$REFERENCE_ALLELE[row], '_', st1$EFFECTIVE_ALLELE[row])])
		} else {
			tmp <- rbind(tmp, shMA[MarkerName == st1$SNP[row]])
			}
	}

for (i in c(1:nrow(tmp))) {
       if (tolower(shGWAS$ra[i]) == tmp$Allele2[i]) {
	       if (tolower(shGWAS$ea[i]) != tmp$Allele1[i]) {
		       print(paste0('Third allele in row #', i))
		       }
	       } else {
		       if (tolower(shGWAS$ra[i]) == tmp$Allele1[i]) {
			        if (tolower(shGWAS$ea[i]) == tmp$Allele2[i]) {
					print(paste0('Flip in row #', i))
					tmp$Allele2[i] <- tolower(shGWAS$ra[i])
					tmp$Allele1[i] <- tolower(shGWAS$ea[i])
					tmp$Effect[i] <- (-1)*tmp$Effect[i]
					tmp$Freq1[i]<- 1 - tmp$Freq1[i]
					dir <- tmp$Direction[i]
					dir <- gsub("\\+", "x",dir)
					dir <- gsub("-", "\\+",dir)
					dir <- gsub("x", "-",dir)
					tmp$Direction[i] <- dir
					} else {
						tmp$Allele1[i] <- tmp$Allele2[i]
						tmp$Allele2[i] <- tolower(shGWAS$ra[i])		  
				  		print(paste0('Third allele in row #', i))
						}
		      } else {
			      print('Error!')
		       }
    }
}
shMA <- tmp
st1 <- mutate(st1, EAF_MA_REPL = shMA$Freq1)
st1 <- mutate(st1, EAF_SE = shMA$FreqSE)
st1 <- mutate(st1, BETA_MA_REPL = shMA$Effect)
st1 <- mutate(st1, SE_MA_REPL = shMA$StdErr)
st1 <- mutate(st1, PVAL_MA_REPL = shMA$'P-value')
st1 <- mutate(st1, Direction = shMA$Direction)
st1 <- mutate(st1, HetISq = shMA$HetISq)
st1 <- mutate(st1, HetChiSq = shMA$HetChiSq)
st1 <- mutate(st1, HetDf = shMA$HetDf)
st1 <- mutate(st1, HetPVal = shMA$HetPVal)
st1 <- mutate(st1, N_total = shMA$n_total)

# Creating an MA-EUR-part of the table

shMAeur <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/03_sh_ma_eur/unification_results/SH_ma_eur_output_done.csv', data.table = F)
ind <- match(st1$SNP, shMAeur$rs_id)
shMAeur <- shMAeur[ind,]
st1 <- mutate(st1, BETA_MA_EUR = shMAeur$beta)
st1 <- mutate(st1, SE_MA_EUR = shMAeur$se)
st1 <- mutate(st1, PVAL_MA_EUR = shMAeur$p)
st1 <- mutate(st1, HWE_MA_EUR = NA)
st1 <- mutate(st1, EAF_MA_EUR = shMAeur$eaf)
st1 <- mutate(st1, N_MA_EUR = shMAeur$n)

# Creating a Back pain-part of the table

st1 <- mutate(st1, BETA_BACK_DISC = backGWAS$beta)
st1 <- mutate(st1, SE_BACK_DISC = backGWAS$se)
st1 <- mutate(st1, PVAL_BACK_DISC = backGWAS$p)
st1 <- mutate(st1, HWE_BACK_DISC = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_BACK_DISC = backGWAS$info)
st1 <- mutate(st1, EAF_BACK_DISC = backGWAS$eaf)
st1 <- mutate(st1, N_BACK_DISC = backGWAS$n)

# Creating a Neck pain-part of the table

neckGWAS <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Neck/Neck_output_done.csv', data.table = F)
ind <- match(st1$SNP, neckGWAS$rs_id)
neckGWAS <- neckGWAS[ind,]
st1 <- mutate(st1, BETA_NECK_DISC = neckGWAS$beta)
st1 <- mutate(st1, SE_NECK_DISC = neckGWAS$se)
st1 <- mutate(st1, PVAL_NECK_DISC = neckGWAS$p)
st1 <- mutate(st1, HWE_NECK_DISC = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_NECK_DISC = neckGWAS$info)
st1 <- mutate(st1, EAF_NECK_DISC = neckGWAS$eaf)
st1 <- mutate(st1, N_NECK_DISC = neckGWAS$n)

# Creating a Headache-part of the table

headGWAS <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Head/Head_output_done.csv', data.table = F)
ind <- match(st1$SNP, headGWAS$rs_id)
headGWAS <- headGWAS[ind,]
st1 <- mutate(st1, BETA_HEAD_DISC = headGWAS$beta)
st1 <- mutate(st1, SE_HEAD_DISC = headGWAS$se)
st1 <- mutate(st1, PVAL_HEAD_DISC = headGWAS$p)
st1 <- mutate(st1, HWE_HEAD_DISC = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_HEAD_DISC = headGWAS$info)
st1 <- mutate(st1, EAF_HEAD_DISC = headGWAS$eaf)
st1 <- mutate(st1, N_HEAD_DISC = headGWAS$n)

# Creating a Hip pain-part of the table

hipGWAS <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Hip/Hip_output_done.csv', data.table = F)
ind <- match(st1$SNP, hipGWAS$rs_id)
hipGWAS <- hipGWAS[ind,]
st1 <- mutate(st1, BETA_HIP_DISC = hipGWAS$beta)
st1 <- mutate(st1, SE_HIP_DISC = hipGWAS$se)
st1 <- mutate(st1, PVAL_HIP_DISC = hipGWAS$p)
st1 <- mutate(st1, HWE_HIP_DISC = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_HIP_DISC = hipGWAS$info)
st1 <- mutate(st1, EAF_HIP_DISC = hipGWAS$eaf)
st1 <- mutate(st1, N_HIP_DISC = hipGWAS$n)

# Creating a Knee pain-part of the table

kneeGWAS <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Knee/Knee_output_done.csv', data.table = F)
ind <- match(st1$SNP, kneeGWAS$rs_id)
kneeGWAS <- kneeGWAS[ind,]
st1 <- mutate(st1, BETA_KNEE_DISC = kneeGWAS$beta)
st1 <- mutate(st1, SE_KNEE_DISC = kneeGWAS$se)
st1 <- mutate(st1, PVAL_KNEE_DISC = kneeGWAS$p)
st1 <- mutate(st1, HWE_KNEE_DISC = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_KNEE_DISC = kneeGWAS$info)
st1 <- mutate(st1, EAF_KNEE_DISC = kneeGWAS$eaf)
st1 <- mutate(st1, N_KNEE_DISC = kneeGWAS$n)

# Creating a Stomach pain-part of the table

stomGWAS <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Stom/Stom_output_done.csv', data.table = F)
ind <- match(st1$SNP, stomGWAS$rs_id)
stomGWAS <- stomGWAS[ind,]
st1 <- mutate(st1, BETA_STOM_DISC = stomGWAS$beta)
st1 <- mutate(st1, SE_STOM_DISC = stomGWAS$se)
st1 <- mutate(st1, PVAL_STOM_DISC = stomGWAS$p)
st1 <- mutate(st1, HWE_STOM_DISC = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_STOM_DISC = stomGWAS$info)
st1 <- mutate(st1, EAF_STOM_DISC = stomGWAS$eaf)
st1 <- mutate(st1, N_STOM_DISC = stomGWAS$n)

fwrite(st1, file = '/mnt/polyomica/projects/bp-sh/data/10_tables/st1_sh_extended.tsv', sep = '\t', dec = ',')









