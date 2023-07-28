invisible(sapply(list.files(pattern = "[.]R$",path="/home/lima/nadya/FREGAT/sumFREGAT_1.2.6.0/R/", full.names = TRUE),source))

#an <- c('nsyn', 'cod','ncod')
an <- 'ncod'
anno <- 'ncod'
if (an == 'cod') anno <- c('nsyn', 'cod')
gene.file <- 'ensembl107.hg38.txt'
cor.path <-'/home/common/DataStorage/UKBB/Project_59345/Matrix/imputed_hg38_white_nonorel_maf-5_machr23_geno02/RData'

tra <- c('SH_disc','BP-SH_disc')
n <- 265000

for (i in 1:2) {

tr <- tra[i]
score.fn <- paste0(tr, '.full.vcf.gz')

tmp <- ACAT(score.file = score.fn, gene.file = gene.file, write.file = paste0('ACAT.', tr,'.',an, '.out'), anno = anno, quiet = T)
tmp <- SKAT(score.file = score.fn,cor.path = cor.path, beta.par = c(1,1), method = "kuonen", rho = TRUE, gene.file = gene.file, gen.var.weights
= 'af', approx = F, write.file = paste0('SKATO.', tr,'.',an, '.out'), anno = anno, p.thr = 0, quiet = T)
tmp <- PCA(score.file = score.fn, cor.path = cor.path, n = n, gene.file = gene.file, anno = anno, gen.var.weights = 'af', approx = F, 
write.file = paste0('PCA.', tr,'.',an,'.out'), quiet = T)

}


