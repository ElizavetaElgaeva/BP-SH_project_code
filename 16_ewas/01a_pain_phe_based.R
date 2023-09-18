#### phenotypes of 4 pains: Neck or shoulder; Back; Hip; Knee pain

library(data.table)

df0 <- fread('ukb52502_6159.txt',h=T,data=F)
df0 <- df0[,1:8]

df <- df0[!is.na(df0$V2),]
exc <- c(-3,8)
inc <- c(3,4,6,7)
v <- which(df$V2 %in% exc)
df <- df[-v,]
for (i in 1:ncol(df)) df[is.na(df[,i]),i] <- 99 

df$phe3 <- df$phe4 <- df$phe6 <- df$phe7 <- 0

for(i in 1:dim(df)[1]) {
for(j in 2:8){
if (df[i,j] == 3) df$phe3 [i] <- 1
if (df[i,j] == 4) df$phe4 [i] <- 1
if (df[i,j] == 6) df$phe6 [i] <- 1
if (df[i,j] == 7) df$phe7 [i] <- 1
}
}

write.table(df, file = 'pain_phe_based.dat', col = T, row = F, qu = F, sep = ' ')

#phe3 <- fread('ukb52502_3404.txt',h=T,data=F)
#phe4 <- fread('ukb52502_3414.txt',h=T,data=F)
#phe6 <- fread('ukb52502_3571.txt',h=T,data=F)
#phe7 <- fread('ukb52502_3773.txt',h=T,data=F)
