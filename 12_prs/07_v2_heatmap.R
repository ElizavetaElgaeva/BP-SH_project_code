# Aim of this script is to create heatmap plot for glm() results for ICD10 and OPCS codes
# representing clusters (h = 1.5) against standardized PRS, test sample, non-relatives 

library(data.table)
library(gplots)

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load table with glm() results
tab <- fread("glm_results_extended_desc_test_nonrelatives.txt", data.table = F)
colnames(tab)

pp <- tab[, c("bp_p", "sh_p", "bp-sh_p")]
thr <- 3.646973e-05
nonsig <- which(pp > thr)

traits <- tab$description

tab <- tab[, c("bp_beta", "sh_beta", "bp-sh_beta")]
tab <- as.matrix(tab)
rownames(tab) <- traits
colnames(tab) <- c("CBP PRS", "SGIT PRS", "CBP UGIT PRS")

tab_p <- round(tab, digits = 2)
tab_p[nonsig] <- NA

# Load cluster data
load("hclust_1.5.RData")

# Select the most representative traits from clusters
table(names(h2) %in% rownames(tab))
which(!(rownames(tab) %in% names(h2)))
#233 234 235 236 237 238
rownames(tab)[c(233, 234, 235, 236, 237, 238)]
table(names(h2) %in% rownames(tab)[-c(233, 234, 235, 236, 237, 238)]) #all true
i <- match(
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
which(rg_p_cut[c6, 1] == min(rg_p_cut[c6, 1]))
i6 <- which(rg_p_cut[c6, 1] == min(rg_p_cut[c6, 1]))
min(rg_p_cut[c6, 2][i6])
traits <- c(traits, names(which(rg_p_cut[c6, 2][i6] == min(rg_p_cut[c6, 2][i6])))) # add trait from the 6th cluster

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
which(rg_p_cut[c11, 1] == min(rg_p_cut[c11, 1]))
i11 <- which(rg_p_cut[c11, 1] == min(rg_p_cut[c11, 1]))
min(rg_p_cut[c11, 2][i11])
traits <- c(traits, names(which(rg_p_cut[c11, 2][i11] == min(rg_p_cut[c11, 2][i11]))))

h2[traits] # self-check, each cluster is represented
rg_m_cut <- rg_m_cut[traits, ]
rg_p_cut <- rg_p_cut[traits, ]

clusters <- c("Pain-related traits and surgery", "Thrombosis", "Injuries", "Respiratory illnesses and smoking", "Anthropometric traits", "Osteoarthritis and other musculoskeletal disorders", "Socio-economic and family status", "Cardio-vascular traits", "Physical activity", "Psychometric traits", "Head pain and medication for pain relief")

rownames(rg_m_cut) <- rownames(rg_p_cut) <- clusters

thr <- 0.05 / (((324*324)-324)/2 + (730 - 322)*2)

rg_m_p <- round(rg_m_cut, digits = 2)
notsig <- which(rg_p_cut > thr)
rg_m_p[notsig] <- NA

pdf("heatmap.pdf", width = 50, height = 50)
heatmap.2(tab, cellnote = tab_p, notecex = 4.0, notecol = "black",
          col = redblue(75), margins = c(25, 120), density.info = "none", trace = "none",
          cexRow = 4.0, srtCol = 45, cexCol = 4.0,
          lhei = c(1, 8), lwid = c(1.2, 4), Colv = FALSE, dendrogram = "row")
dev.off()

save(traits, clusters, rg_m_cut, rg_p_cut, file = "heatmap_clusters_1.8.RData")
