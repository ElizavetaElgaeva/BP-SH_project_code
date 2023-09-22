invisible(sapply(list.files(pattern = "[.]R$",path="/home/lima/nadya/FREGAT/sumFREGAT_1.2.6_collaps/R/", full.names = TRUE), source))
tr <- 'SGIT'
mac <- 'mac10'
meth <- c('SKATO','PCA')
ann <- c('lof','nsyn','cod','all')
score.file <- paste(tr,mac,'full.vcf.gz',sep='.')
gene.file <- '../ensembl107.hg38.txt'
cor.path <- '/home/common/DataStorage/UKBB/Project_59345/Matrix/exome_white_qc_nonrel_mac3_geno02/RData/'
gene <- 'all'

for (an in ann) {
if (an == 'lof'){
anno <- 'lof'
set <- 'set1'
}
if (an == 'nsyn'){
anno <- c('lof','nsyn')
set <- 'set2'
}
if (an == 'cod'){
anno <- c('lof','nsyn','cod')
set <- 'set3'
}
if (an == 'all'){
anno <- c('lof','nsyn','cod','ncod')
set <- 'set4'
}

collaps.fn <- paste0('/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/',tr,'/',tr,'_self-rep_eur_exomes_4_pains_',set,'.',mac,'_collaps.fastGWA')
col <- read.table(collaps.fn, h = T)
colnames(col) <- toupper(colnames(col))
colnames(col)[colnames(col) == 'B'] <- 'BETA'
colnames(col)[colnames(col) == 'EAF'] <- 'AF1'
collaps.fn <- paste(tr,an,'tmp',sep='.')
write.table(col, file = collaps.fn, qu = F, row = F)

for (m in meth) {
flout <- paste(tr,m,an,mac,sep='.')

if (m == 'SKATO') {
tmp <- SKAT(score.file = score.file, genes = gene, gene.file = gene.file, cor.path = cor.path, anno = anno, collaps.fn = collaps.fn, write.file = flout, quiet = T)
} else {
tmp <- PCA(score.file = score.file, genes = gene, gene.file = gene.file, cor.path = cor.path, anno = anno, collaps.fn = collaps.fn, n=179000,write.file = flout, quiet = T)
}
}
}

