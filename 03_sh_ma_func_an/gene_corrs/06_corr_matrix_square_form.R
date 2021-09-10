library(data.table)
library(tidyr)

rg_results <- fread('../../../../data/03_sh_ma_func_an/gene_corrs/sq_matrix_long_format.txt', data.table=F)
se_res <- rg_results[,c(1,2,6)]
rg_results <- rg_results[,c(1,2,4)]

rg_results_double <- rg_results[,c(2,1,3)]
se_res_double <- se_res[,c(2,1,3)]

colnames(rg_results_double)<-colnames(rg_results)
colnames(se_res_double)<-colnames(se_res)


rg_results <- rbind(rg_results, rg_results_double)
se_res <- rbind(se_res, se_res_double)

wide_rg <- spread(rg_results, gwas_id_2, rg)
wide_se <- spread(se_res, gwas_id_2, se)

'Number of NA should be:'
nrow(wide_rg)
'Actual number:'
length(wide_rg[is.na(wide_rg)])

wide_rg[is.na(wide_rg)]<-1
wide_se[is.na(wide_se)]<-0

rownames(wide_rg)<-wide_rg$gwas_id_1
rownames(wide_se)<-wide_se$gwas_id_1
wide_rg<-wide_rg[,-1]
wide_se<-wide_se[,-1]
write.csv(wide_rg,'../../../../data/03_sh_ma_func_an/gene_corrs/gene_corr_matrix.csv',row.names=T)
write.csv(wide_se,'../../../../data/03_sh_ma_func_an/gene_corrs/gene_corr_se_matrix.csv',row.names=T)

