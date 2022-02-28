# im of this script is to perform preliminary clumping in adjusted BP GWAS (Discovery)
# to check whether there is any signal

library(data.table)

setwd("/mnt/polyomica/projects/bp-sh/src/elgaeva_src/06_clumping")
source("./01_sh/f_clumping.R")

dt <- fread("/mnt/polyomica/projects/bp-sh/data/13_bp_adj/bp_adj_discovery/BP_adj_discovery.txt", data.table = F)
colnames(dt)
thr <- 5e-8/3
res <- function_for_shlop_29_03_2020(dt, p_value = "p", pos = "pos", snp = "SNP", thr = thr)
res

