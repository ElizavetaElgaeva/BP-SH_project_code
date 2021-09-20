# Aim of this script is to perform hierarchical clustering of the traits 
# statistically significantly genetically correlated with either SH or BP-SH
# with |rg| > 0.25

library(data.table)

setwd("/mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/gene_corrs")

rg_m <- fread("./gene_corr_matrix.csv", data.table = F)
rownames(rg_m) <- rg_m[1, ]
colnames(rg_m) <- rg_m[1, ]
rg_m <- rg_m[-1, -1] # delete row and column with gwas_id
rg_m_cut <- rg_m[-c(323:324), -c(323:324)] # keep only traits correlated with SH and BP-SH

desc <- fread("/mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi/20210901_DESCRIPTORS.csv", data.table = F)
ind <- match(colnames(rg_m_cut), desc$gwas_id)
desc <- desc[ind, ]

rownames(rg_m_cut) <- colnames(rg_m_cut) <- desc$trait_name

dd <- (1 - abs(rg_m_cut))

d <- as.dist(dd)
h1 <- hclust(d, method = "ward.D2")

pdf("hclust.pdf", width = 70, height = 70)
plot(h1, hang = -1)
abline(h = 1.8, col = "blue")
abline(h = 1.7, col = "red")
abline(h = 1.5, col = "green")
dev.off()

library(corrplot)

desc <- fread("/mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi/20210901_DESCRIPTORS.csv", data.table = F)
ind2 <- match(colnames(rg_m), desc$gwas_id)
desc2 <- desc[ind2, ]

#rownames(rg_m) <- colnames(rg_m) <- desc2$trait_name
#rg_m[rg_m > 1] <- 1
#rg_m[rg_m < -1] <- -1
rg2plot <- rg_m_cut[h1$order, h1$order]

pdf("heatmap_hclust.pdf", width = 70, height = 70)
corrplot(as.matrix(rg2plot), type = 'full', is.corr = F)
dev.off()

h2 <- cutree(tree = h1, h = 1.8)
save(h2, file = "hclust_1.8.RData")
