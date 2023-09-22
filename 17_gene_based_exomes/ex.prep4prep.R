#### prep4prep
library(data.table)
tr <- 'SGIT'
mac <- 10
dir <- '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/'
fl <- paste0(dir,tr,'/',tr,'_self-rep_eur_exomes_4_pains.fastGWA') #gwas
df <- fread(fl, h = T, data.table = F)
df$MAC <- round(df$eaf * df$N_mean *2)
df <- df[pmin(df$eaf,1-df$eaf)  <= 0.01 & df$MAC > mac,] #filters: MAC > mac & maf <=0.01 

ref <- get(load('/home/common/DataStorage/UKBB/Project_59345/Matrix/exome_white_qc_nonrel_mac3_geno02/RData/reference.RData'))
ref <- ref[ref$AF <= 0.01,]

df$ID <- paste(df$SNP,df$A1,df$A2,sep=':')  ## name as in reference file

df <- df[df$ID %in% ref$ID,] ##markers in reference file

vep <- fread('~/Zor/BACK_PAIN/Exome_VEP_annotation/exome.anno.single.txt', h = T, data.table = F)
df$ind <- paste(df$chr,df$pos,df$A1,df$A2,sep=':')  ## name as in vep

v <- match(vep$snp,df$ind)
vep$ID <- df$ID[v]  ##name as in reference file
vep$ALT <- df$A1[v]
vep$P <- df$p[v]
vep$BETA <- df$b[v]
vep$AF1 <- df$eaf[v]

vep <- na.omit(vep)


invisible(sapply(list.files(pattern = "[.]R$", path="/home/lima/nadya/FREGAT/sumFREGAT_1.2.6.0/R/", full.names = TRUE), source))
ref <- get(load('/home/common/DataStorage/UKBB/Project_59345/Matrix/exome_white_qc_nonrel_mac3_geno02/RData/reference.RData'))
input.data <- data.frame(CHROM = vep$chr, POS = vep$pos, ID = vep$ID, EA = vep$ALT, P = vep$P, BETA = vep$BETA, EAF = vep$AF1, GENE = vep$gene, ANNO = vep$newANNO)
prep.score.files(input.data, ref = ref, output = paste0(tr, '.mac',mac,'.full'))


