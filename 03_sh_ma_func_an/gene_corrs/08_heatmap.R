# Aim of this script is to create heatmap plot for traits representing clusters (h = 1.8),
# obtained at previous step

library(data.table)
library(gplots)

setwd("/mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/gene_corrs")

# Load rg matrix
rg_m <- fread("./gene_corr_matrix.csv", data.table = F)
rg_m <- as.matrix(rg_m)
rownames(rg_m) <- rg_m[1, ]
colnames(rg_m) <- rg_m[1, ]
rg_m <- rg_m[-1, -1] # delete row and column with gwas_id
rg_m_cut <- rg_m[-c(323:324), c(323:324)] # keep only traits correlated with SH and BP-SH

# load p-value matrix
rg_p <- fread("./gene_corr_p_matrix.csv", data.table = F)
rg_p <- as.matrix(rg_p)
rownames(rg_p) <- rg_p[1, ]
colnames(rg_p) <- rg_p[1, ]
rg_p <- rg_p[-1, -1] # delete row and column with gwas_id
rg_p_cut <- rg_p[-c(323:324), c(323:324)] # keep only traits correlated with SH and BP-SH

# Load cluster data
load("hclust_1.8.RData")

# Load descriptors
desc <- fread("/mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi/20210901_DESCRIPTORS.csv", data.table = F)

# Rename
ind1 <- match(rownames(rg_m_cut), desc$gwas_id)
rownames(rg_m_cut) <- rownames(rg_p_cut) <- desc$trait_name[ind1]
colnames(rg_m_cut) <- colnames(rg_p_cut) <- c("SGIT", "CBP UGIT")

# Select the most representative traits from clusters
table(names(h2) == rownames(rg_m_cut)) # all true
table(h2)
#h2
# 1  2  3  4  5  6  7  8  9 10 11
#70  8 24 31 37 14 30 27 48 24  9
h2 <- sort(h2, decreasing = FALSE)
rg_m_cut <- rg_m_cut[names(h2), ]
rg_p_cut <- rg_p_cut[names(h2), ]

c1 <- which(h2 == 1)
min(rg_p_cut[c1, ])
i1 <- which(rg_p_cut == min(rg_p_cut[c1, ]))
summary(rg_p_cut[i1, 1]) # all 0
summary(rg_p_cut[i1, 2])
traits <- c(names(which(rg_p_cut[i1, 2] == min(rg_p_cut[i1, 2])))) # trait from the 1st cluster

c2 <- which(h2 == 2)
min(rg_p_cut[c2, ])
which(rg_p_cut == min(rg_p_cut[c2, ]))
i2 <- which(rg_p_cut == min(rg_p_cut[c2, ]))
traits <- c(traits, names(h2)[i2]) # add trait from the 2nd cluster

c3 <- which(h2 == 3)
min(rg_p_cut[c3, ])
which(rg_p_cut == min(rg_p_cut[c3, ]))
i3 <- which(rg_p_cut == min(rg_p_cut[c3, ]))
traits <- c(traits, names(h2)[i3]) # add trait from the 3d cluster

c4 <- which(h2 == 4)
min(rg_p_cut[c4, ])
which(rg_p_cut == min(rg_p_cut[c4, ]))
i4 <- which(rg_p_cut == min(rg_p_cut[c4, ]))
traits <- c(traits, names(h2)[i4]) # add trait from the 4th cluster

c5 <- which(h2 == 5)
min(rg_p_cut[c5, ])
which(rg_p_cut == min(rg_p_cut[c5, ]))
i5 <- which(rg_p_cut == min(rg_p_cut[c5, ]))
traits <- c(traits, names(h2)[i5]) # add trait from the 5th cluster

c6 <- which(h2 == 6)
min(rg_p_cut[c6, ])
which(rg_p_cut == min(rg_p_cut[c6, ]))
i6 <- which(rg_p_cut[, 2] == min(rg_p_cut[c6, 2]))
traits <- c(traits, names(h2)[i6]) # add trait from the 6th cluster

c7 <- which(h2 == 7)
min(rg_p_cut[c7, ])
which(rg_p_cut == min(rg_p_cut[c7, ]))
i7 <- which(rg_p_cut == min(rg_p_cut[c7, ]))
traits <- c(traits, names(h2)[i7]) # add trait from the 7th cluster

c8 <- which(h2 == 8)
min(rg_p_cut[c8, ])
which(rg_p_cut == min(rg_p_cut[c8, ]))
i8 <- which(rg_p_cut == min(rg_p_cut[c8, ]))
traits <- c(traits, names(h2)[i8]) # add trait from the 8th cluster

c9 <- which(h2 == 9)
min(rg_p_cut[c9, ])
which(rg_p_cut == min(rg_p_cut[c9, ]))
i9 <- which(rg_p_cut == min(rg_p_cut[c9, ]))
traits <- c(traits, names(h2)[i9]) # add trait from the 9th cluster

c10 <- which(h2 == 10)
min(rg_p_cut[c10, ])
which(rg_p_cut == min(rg_p_cut[c10, ]))
i10 <- which(rg_p_cut == min(rg_p_cut[c10, ]))
traits <- c(traits, names(h2)[i10]) # add trait from the 10th cluster

c11 <- which(h2 == 11)
min(rg_p_cut[c11, ])
which(rg_p_cut == min(rg_p_cut[c11, ]))
i11 <- which(rg_p_cut[, 2] == min(rg_p_cut[c11, 2]))
traits <- c(traits, names(h2)[i11]) # add trait from the 11th cluster

h2[traits] # self-check, each cluster is represented
rg_m_cut <- rg_m_cut[traits, ]
rg_p_cut <- rg_p_cut[traits, ]

clusters <- c("Pain-related traits and surgery", "Thrombosis", "Injuries", "Respiratory illnesses and smoking", "Anthropometric traits", "Osteoarthritis and other musculoskeletal disorders", "Socio-economic and family status", "Cardio-vascular traits", "Physical activity", "Psychometric traits", "Head pain and medication for pain relief")

rownames(rg_m_cut) <- rownames(rg_p_cut) <- clusters

thr <- 0.05 / (((324*324)-324)/2 + (730 - 322)*2)

rg_m_p <- round(rg_m_cut, digits = 2)
notsig <- which(rg_p_cut > thr)
rg_m_p[notsig] <- NA

pdf("heatmap_clust_1.8.pdf", width = 30, height = 30)
heatmap.2(rg_m_cut, cellnote = rg_m_p, notecex = 4.0, notecol = "black", na.color = "grey",
	  col = redblue(75), margins = c(25, 80), density.info = "none", trace = "none",
          cexRow = 4.0, srtCol = 45, cexCol = 4.0,
	  lhei = c(1, 8), lwid = c(1.2, 4), Colv = FALSE, dendrogram = "row")
dev.off()

save(traits, clusters, rg_m_cut, rg_p_cut, file = "heatmap_clusters_1.8.RData")
