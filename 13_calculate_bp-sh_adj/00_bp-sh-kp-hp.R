# Aim of this script is to calculate the linear combination coefficients for UGIT CBP 
# adjusted for chronic headache and knee pain 

#functions

source("../Dropbox/paul_aging/leve_one_out_gip/00_core_functions.R")

ugct_coefficients_estimation=function(gcov,phe,y,ind){
  #gcov - genetic covariance matrix
  #phe - phenotypic correlashions
  #y - index of the trait that will be adjusted
  #ind - indexes of traits that will be used for adjustment
  
  gcov=gcov[c(y,ind),c(y,ind)]
  phe=phe[c(y,ind),c(y,ind)]
  
  b=.b(response = 1,pred = 2:nrow(gcov),S=gcov)[,1]
  #.t(response = 1,pred = 2:5,S=gcov,N=1e5)
  ugct=c(1,rep(0,length(b)))-c(0,b)
  names(ugct)=colnames(gcov)
  v1=phen(a=ugct,phem=phe)
  ugct=ugct/sqrt(v1)
  
  exp_h2=H2(a = ugct,covm = gcov,phem = phe)
  exp_gcor_with_orig=cor_gi_alfa(a=ugct,i=1,covm=gcov)
  output=list(coeffs=ugct,expected_h2=exp_h2,expected_gcor_with_orig=exp_gcor_with_orig)
  return(output)
}


###########
setwd("../Dropbox/BP-SH/results/BP-SH-KP/")

gcov=read.table("gene_cov_matrix.txt",header=T,row.names = 1)
pcov=read.table("pheno_corr_matrix.txt",header=T,row.names = 1)
colnames(gcov)=rownames(gcov)=colnames(pcov)=rownames(pcov)=
  c("Hip","Back","Neck", "Knee", "Head", "Stomach", "SGC", "UGC of CBP")


y=match("UGC of CBP",colnames(gcov))

ind=match(c("Hip","Knee", "Head", "Stomach"),colnames(gcov))
out=ugct_coefficients_estimation(gcov=as.matrix(gcov),phe=as.matrix(pcov),y=y,ind=ind)
out

ind=match(c("Knee","Head"),colnames(gcov))
out=ugct_coefficients_estimation(gcov=as.matrix(gcov),phe=as.matrix(pcov),y=y,ind=ind)
out

H2(a=out$coeffs,covm = as.matrix(gcov)[names(out$coeffs),names(out$coeffs)],phem = as.matrix(pcov)[names(out$coeffs),names(out$coeffs)])

cor_gi_alfa(a=out$coeffs,i=3,covm=as.matrix(gcov)[names(out$coeffs),names(out$coeffs)])
out

