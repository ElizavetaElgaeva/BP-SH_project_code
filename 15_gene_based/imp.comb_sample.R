## Combined three metods by ACATO and combined with collaps by FISHER for every gene set
library(data.table)

trait <- c('SH','BP-SH')
set <- c('nsyn','cod','ncod')
#set <- 'ncod'
smpl <- c('disk','rep','ma')

for (t in trait) {
for (s in set) {
fl <- paste(t,s,'comb.csv',sep='.')
df <- read.csv2(fl,h=T)
df <- df[,c(1:5,9)]
for (i in 2:3) {
m <- smpl[i]
if(i == 2) fl <- paste0('./SH_rep_on_replication_sample/ALLgenes/',t,'.',s,'.comb.csv')
if(i == 3) fl <- paste0('./SH_rep_on_ma_sample/ALLgenes/',t,'.',s,'.comb.csv')
df1 <- read.csv2(fl,h=T)
df <- cbind(df,df1$markers)
colnames(df)[dim(df)[2]] <- paste0('markers.',m)
df <- cbind(df,df1$COMB)
colnames(df)[dim(df)[2]] <- paste0('COMB.',m)
}

for (i in 1:dim(df)[1]){
df$p.min[i] <- min(c(df[i,6],df[i,8],df[i,10]))
}

## Remove  NA
#df <- df[!is.na(df$p.min),]

fl <- paste(t,s,'comb.sample.csv',sep='.')
write.csv2(df,fl,row.names = F,quote=F)
fl <- paste(t,s,'comb.sample.sig.csv',sep='.')
sig <- na.omit(df[df$p.min <= 2.5E-06,])
write.csv2(sig,fl,row.names = F,quote=F)

}
}
