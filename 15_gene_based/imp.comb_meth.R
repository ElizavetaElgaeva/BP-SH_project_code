## Combined three metods by ACATO and combined with collaps by FISHER for every gene set
library(data.table)
invisible(sapply(list.files(pattern = "[.]R$", path="/home/lima/nadya/FREGAT/sumFREGAT_1.2.6.0/R/", full.names = TRUE), source))

df0 <- fread('ensembl107.hg38.txt',h=F,data=F)
df0 <- df0[,c(1,3,5,6)]
colnames(df0) <- c('gene','chr','start','stop')

trait <- c('SH','BP-SH')
set <- c('nsyn','cod','ncod')
meth <- c('ACAT','SKATO','PCA')

for (t in trait) {
for (s in set) {
df <- df0
for (m in meth) {
fl <- paste(m,t,s,'out',sep='.')
df1 <- fread(fl,h=T,data=F)
df <- cbind(df,df1$pvalue)
colnames(df)[dim(df)[2]] <- m
}

df$markers <- df1$markers
df$filtered.markers <- df1$filtered.markers

for (i in 1:dim(df)[1]){
df$COMB[i] <- ACATO(c(df$ACAT[i],df$SKATO[i],df$PCA[i]))
}

## Remove  NA
df <- df[!is.na(df$COMB),]

fl <- paste(t,s,'comb.csv',sep='.')
write.csv2(df[,c(1:4,8,5:7,10)],fl,row.names = F,quote=F)
fl <- paste(t,s,'comb.sig.csv',sep='.')
sig <- na.omit(df[df$COMB <= 2.5E-06,])
write.csv2(sig[,c(1:4,8,5:7,10)],fl,row.names = F,quote=F)

}
}
