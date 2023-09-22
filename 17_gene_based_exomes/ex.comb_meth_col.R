## Combined three metods by ACATO and combined with collaps by FISHER for every gene set
library(data.table)
invisible(sapply(list.files(pattern = "[.]R$", path="/home/lima/nadya/FREGAT/sumFREGAT_1.2.6.0/R/", full.names = TRUE), source))

t <- 'SGIT'

df0 <- fread('../ensembl107.hg38.txt',h=F,data=F)
df0 <- df0[,c(1,3,5,6)]
colnames(df0) <- c('gene','chr','start','stop')
dir <- '/home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CBP/'

set <- c('lof','nsyn','cod','all')
mac <- 'mac10'
meth <- c('SKATO','PCA')

flsig <- paste(t,'exome.sig.csv',sep='.')
sig.comb <- c()

for (s in set) {
for (mc in mac) {
df <- df0
for (m in meth) {
fl <- paste(t,m,s,mc,sep='.')
df1 <- fread(fl,h=T,data=F)
df <- cbind(df,df1$pvalue)
colnames(df)[dim(df)[2]] <- m
}

df$markers <- df1$markers
df$filtered.markers <- df1$filtered.markers

for (i in 1:dim(df)[1]){
df$COMB[i] <- ACATO(c(df$SKATO[i],df$PCA[i]))
}

## Remove  NA
df <- df[!is.na(df$COMB),]
df$set <- s
fl <- paste(t,s,mc,'csv',sep='.')
#write.csv2(df[,c(1:4,7,8,5,6,9)],fl,row.names = F,quote=F)

sig <- na.omit(df[df$COMB <= 2.5E-05,])
sig.comb <- rbind(sig.comb,sig)
}
}
write.csv2(sig.comb[,c(1:4,7,8,5,6,9,10)],flsig,row.names = F,quote=F)

