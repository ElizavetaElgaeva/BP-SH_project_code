invisible(sapply(list.files(pattern = "[.]R$",path="/home/lima/nadya/FREGAT/sumFREGAT_1.2.6.0/R/", full.names = TRUE),source))
anno <- 'ncod'
#anno <- 'nsyn'
#anno <- c('nsyn', 'cod')
an <- 'ncod'
gene.file <- 'hg38.107'
cor.path <-'/home/common/DataStorage/UKBB/Project_59345/Matrix/imputed_hg38_white_nonorel_maf-5_machr23_geno02/RData'

tr0 <- c('SH', 'BP-SH')

n <- 265000

for (i in 1:1) {

tr <- tr0[i]
score.file <- paste0(tr, '_disc.full.vcf.gz')

tmp <- ACAT(score.file = score.file, gene.file = gene.file, write.file = paste0('ACAT.', tr,'.',an, '.out'), anno = anno, quiet = T)
#tmp <- SKAT(score.file = score.file,cor.path = cor.path, beta.par = c(1,1), method = "kuonen", rho = TRUE, gene.file = gene.file, gen.var.weights
#= 'af', approx = F, write.file = paste0('SKATO.', tr,'.',an, '.out'), anno = anno, p.thr = 0, quiet = T)
#tmp <- PCA(score.file = score.file, cor.path = cor.path, n = n, gene.file = gene.file, anno = anno, gen.var.weights = 'af', approx = F, write.file
#= paste0('PCA.', tr,'.',an,'.out'), quiet = T)

}


