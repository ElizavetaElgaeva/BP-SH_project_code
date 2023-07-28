#### prep4prep
library(data.table)
tr <- 'BP-SH_disc'
fl <- paste0(tr,'_output_done.csv') #gwas
df0 <- fread(fl, h = T, data.table = F)
df <- df0[,c(3,7:13)]

an0 <- get(load('gs.RDa')) ##hg38
ref <- get(load('/home/common/DataStorage/UKBB/Project_59345/Matrix/imputed_hg38_white_nonorel_maf-5_machr23_geno02/RData/reference.RData')) ##hg38
## 2 new IDs
v <- which(!grepl('rs', an0$newID))
an <- an0[-v, ] ##only with rs
tmp <- strsplit(an$newID, ':')
tmp <- unlist(tmp)
rs <- tmp[grepl('rs', tmp, fixed = T)]
an$rs_id <- rs ##rs name

dat <- merge(an,df, by = "rs_id")
#dim(dat)
#[1] 5588236      16


input.data <- data.frame(CHROM = dat$chr, POS = dat$POS_GB, ID = dat$newID, EA = dat$ea, P = dat$p, BETA = dat$beta, EAF = dat$eaf, GENE = dat$gene, ANNO = dat$anno)
write.table(input.data, file = paste0(tr,'.full'),row = F, col = T, qu = F, sep = ' ')

invisible(sapply(list.files(pattern = "[.]R$", path="/home/lima/nadya/FREGAT/sumFREGAT_1.2.6.0/R/", full.names = TRUE), source))
prep.score.files(input.data, ref = ref, output = paste0(tr,'.full'))

