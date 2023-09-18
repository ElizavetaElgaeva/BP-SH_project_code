## create pheno_chron

library(data.table)
df0 <- fread('pain_phe_based.dat',,h=T,data=F)
df <- df0[,c(1,9:12)]

#phe0 <- fread('ukb52502_3404.txt',h=T,data=F) ##phe3
phe0 <- fread('ukb52502_3414.txt',h=T,data=F) ## phe6
#phe0 <- fread('ukb52502_3773.txt',h=T,data=F) ## phe7
#phe0 <- fread('ukb52502_3571.txt',h=T,data=F) ##phe4 (back pain)

## names(df): ID phe7 phe6 phe4 phe3
df <- df[,c(1,3)]####
colnames(df)[2] <- 'phe'
phe <- phe0[,1:2]
phe <- phe[!is.na(phe$V2),]

joint <- merge(df,phe,by = 'ID', all.x = T)
joint[is.na(joint$V2),3] <- 99

joint$chron <- 99
for(i in 1:dim(joint)[1]) {
if(joint$phe [i] == 0) joint$chron [i] <- 0
if (joint$phe [i] == 1 & joint$V2 [i] == 0) joint$chron [i] <- 0
if (joint$phe [i] == 1 & joint$V2 [i] == 1) joint$chron [i] <- 1
}
write.table(joint, file = 'pain_phe6_chron.dat', col = T, row = F, qu = F, sep = ' ')
