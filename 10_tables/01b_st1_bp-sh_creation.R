# Aim of this script is to create a full table representing all associated with BP-SH loci
# and provide full information for inference that they are replicated

library(data.table)
library(dplyr)

setwd('/mnt/polyomica/projects/bp-sh/data')

# Creating an SNP-information-part of a table
backGWAS <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/00_original_traits/unification_results/Back/Back_output_done.csv', data.table = F)
bp_shGWAS <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/02_bp-sh/unification_results/BP-SH_disc_output_done.csv', data.table = F)
cojoFile <- fread(input = '/mnt/polyomica/projects/bp-sh/data/07_cojo/02_bp-sh/bp-sh_disc_all_chr_2.5e-8.jma.cojo', data.table = F)
st1 <- cojoFile$SNP
st1 <- as.data.frame(st1, stringsAsFactors=FALSE)
colnames(st1) <- c('SNP')
ind <- match(st1$SNP, bp_shGWAS$rs_id)
bp_shGWAS <- bp_shGWAS[ind,]
st1 <- mutate(st1, CHR = cojoFile$Chr)
st1 <- mutate(st1, POS = cojoFile$bp)
st1 <- mutate(st1, REFERENCE_ALLELE = bp_shGWAS$ra)
st1 <- mutate(st1, EFFECTIVE_ALLELE = bp_shGWAS$ea)

# Creating a discovery-part of the table

ind <- match(st1$SNP, backGWAS$rs_id)
backGWAS <- backGWAS[ind,]
intercept <- 0.987876106164573 # LD Score intercept for BP-SH Discovery calculated using --h2 option through the data base
st1 <- mutate(st1, MAF_DISC = pmin(bp_shGWAS$eaf, (1 - bp_shGWAS$eaf)))
st1 <- mutate(st1, EAF_DISC = bp_shGWAS$eaf)
st1 <- mutate(st1, IMPUTATION_QUALITY_DISC = backGWAS$info)
st1 <- mutate(st1, BETA_DISC = bp_shGWAS$beta)
st1 <- mutate(st1, SE_DISC = bp_shGWAS$se)
st1 <- mutate(st1, PVAL_DISC = bp_shGWAS$p)
st1 <- mutate(st1, PVAL_GC_DISC = pchisq(((st1$BETA_DISC / st1$SE_DISC)^2 / intercept), df = 1, lower.tail = F))
st1 <- mutate(st1, HWE_DISC = NA)
st1 <- mutate(st1, N_DISC = bp_shGWAS$n)

# Creating a replication-EUR-part of the table

bp_shEUR <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/EA/02_bp-sh/unification_results/BP-SH_ea_repl_output_done.csv', data.table = F)
ind <- match(st1$SNP, bp_shEUR$rs_id)
bp_shEUR <- bp_shEUR[ind,]
backEUR <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/EA/00_original_traits/unification_results/Back/Back_output_done.csv', data.table = F)
ind <- match(st1$SNP, backEUR$rs_id)
backEUR <- backEUR[ind,]
st1 <- mutate(st1, BETA_EUR = bp_shEUR$beta)
st1 <- mutate(st1, SE_EUR = bp_shEUR$se)
st1 <- mutate(st1, PVAL_EUR = bp_shEUR$p)
st1 <- mutate(st1, HWE_EUR = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_EUR = backEUR$info)
st1 <- mutate(st1, EAF_EUR = bp_shEUR$eaf)
st1 <- mutate(st1, N_EUR = bp_shEUR$n)

# Creating a replication-AFR-part of the table

bp_shAFR <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/AA/02_bp-sh/unification_results/BP-SH_aa_repl_output_done.csv', data.table = F)
ind <- match(st1$SNP, bp_shAFR$rs_id)
bp_shAFR <- bp_shAFR[ind,]
backAFR <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/AA/00_original_traits/unification_results/Back/Back_output_done.csv', data.table = F)
ind <- match(st1$SNP, backAFR$rs_id)
backAFR <- backAFR[ind,]
st1 <- mutate(st1, BETA_AFR = bp_shAFR$beta)
st1 <- mutate(st1, SE_AFR = bp_shAFR$se)
st1 <- mutate(st1, PVAL_AFR = bp_shAFR$p)
st1 <- mutate(st1, HWE_AFR = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_AFR = backAFR$info)
st1 <- mutate(st1, EAF_AFR = bp_shAFR$eaf)
st1 <- mutate(st1, N_AFR = bp_shAFR$n)

# Creating a replication-ASI-part of the table

bp_shASI <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/SA/02_bp-sh/unification_results/BP-SH_sa_repl_output_done.csv', data.table = F)
ind <- match(st1$SNP, bp_shASI$rs_id)
bp_shASI <- bp_shASI[ind,]
backASI <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/SA/00_original_traits/unification_results/Back/Back_output_done.csv', data.table = F)
ind <- match(st1$SNP, backASI$rs_id)
backASI <- backASI[ind,]
st1 <- mutate(st1, BETA_ASI = bp_shASI$beta)
st1 <- mutate(st1, SE_ASI = bp_shASI$se)
st1 <- mutate(st1, PVAL_ASI = bp_shASI$p)
st1 <- mutate(st1, HWE_ASI = NA)
st1 <- mutate(st1, IMPUTATION_QUALITY_ASI = backASI$info)
st1 <- mutate(st1, EAF_ASI = bp_shASI$eaf)
st1 <- mutate(st1, N_ASI = bp_shASI$n)

# Creating an MA-replication-part of the table

bp_shMA <- fread(input = '/mnt/polyomica/projects/bp-sh/data/05_bp-sh_ma/02_repl_ma/02_ma_METAL/BP-SH_ma_repl_metal.txt')
tmp <- NULL
for (row in c(1:nrow(st1))) {
	if (nrow(bp_shMA[MarkerName == st1$SNP[row]]) == 0) {
		tmp <- rbind(tmp, bp_shMA[MarkerName == paste0(st1$CHR[row], ':', st1$POS[row], '_', st1$REFERENCE_ALLELE[row], '_', st1$EFFECTIVE_ALLELE[row])])
		} else {
			tmp <- rbind(tmp, bp_shMA[MarkerName == st1$SNP[row]])
			}
	}

for (i in c(1:nrow(tmp))) {
       if (tolower(bp_shGWAS$ra[i]) == tmp$Allele2[i]) {
	       if (tolower(bp_shGWAS$ea[i]) != tmp$Allele1[i]) {
		       print(paste0('Third allele in row #', i))
		       }
	       } else {
		       if (tolower(bp_shGWAS$ra[i]) == tmp$Allele1[i]) {
			        if (tolower(bp_shGWAS$ea[i]) == tmp$Allele2[i]) {
					print(paste0('Flip in row #', i))
					tmp$Allele2[i] <- tolower(bp_shGWAS$ra[i])
					tmp$Allele1[i] <- tolower(bp_shGWAS$ea[i])
					tmp$Effect[i] <- (-1)*tmp$Effect[i]
					tmp$Freq1[i]<- 1 - tmp$Freq1[i]
					dir <- tmp$Direction[i]
					dir <- gsub("\\+", "x",dir)
					dir <- gsub("-", "\\+",dir)
					dir <- gsub("x", "-",dir)
					tmp$Direction[i] <- dir
					} else {
						tmp$Allele1[i] <- tmp$Allele2[i]
						tmp$Allele2[i] <- tolower(bp_shGWAS$ra[i])		  
				  		print(paste0('Third allele in row #', i))
						}
		      } else {
			      print('Error!')
		       }
    }
}
bp_shMA <- tmp
st1 <- mutate(st1, EAF_MA_REPL = bp_shMA$Freq1)
st1 <- mutate(st1, EAF_SE = bp_shMA$FreqSE)
st1 <- mutate(st1, BETA_MA_REPL = bp_shMA$Effect)
st1 <- mutate(st1, SE_MA_REPL = bp_shMA$StdErr)
st1 <- mutate(st1, PVAL_MA_REPL = bp_shMA$'P-value')
st1 <- mutate(st1, Direction = bp_shMA$Direction)
st1 <- mutate(st1, HetISq = bp_shMA$HetISq)
st1 <- mutate(st1, HetChiSq = bp_shMA$HetChiSq)
st1 <- mutate(st1, HetDf = bp_shMA$HetDf)
st1 <- mutate(st1, HetPVal = bp_shMA$HetPVal)
st1 <- mutate(st1, N_total = bp_shMA$n_total)

# Creating an MA-EUR-part of the table

bp_shMAeur <- fread(input = '/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/05_bp-sh_ma_eur/unification_results/BP-SH_ma_eur_output_done.csv', data.table = F)
ind <- match(st1$SNP, bp_shMAeur$rs_id)
bp_shMAeur <- bp_shMAeur[ind,]
st1 <- mutate(st1, BETA_MA_EUR = bp_shMAeur$beta)
st1 <- mutate(st1, SE_MA_EUR = bp_shMAeur$se)
st1 <- mutate(st1, PVAL_MA_EUR = bp_shMAeur$p)
st1 <- mutate(st1, HWE_MA_EUR = NA)
st1 <- mutate(st1, EAF_MA_EUR = bp_shMAeur$eaf)
st1 <- mutate(st1, N_MA_EUR = bp_shMAeur$n)

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

fwrite(st1, file = '/mnt/polyomica/projects/bp-sh/data/10_tables/st1_bp-sh_extended.tsv', sep = '\t', dec = ',')



