# Aim of this script is to create heatmap plot for glm() results for ICD10 and OPCS codes
# representing clusters (h = 1.25) against standardized PRS, test sample, non-relatives
# Note: all codes combined to level 2; ICD10 codes from chapters I-XVII only; OPCS codes excluding X, Y, Z chapters

library(data.table)
library(gplots)

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load table with glm() results
tab <- fread("glm_results_level_2_preselected_codes_extended_desc_test_nonrelatives.txt", data.table = F)
colnames(tab)

pp <- tab[, c("bp_p", "sh_p", "bp-sh_p")]
pp <- as.matrix(pp)
colnames(pp) <- c("bp_p", "sh_p", "bp-sh_p")
thr <- 0.05/(3*(165 + 132 + 6)) # 165 is a number of icd10 codes, 132 is a number of opcs codes, 6 is a number of pain traits
nonsig <- which(pp > thr)

traits <- tab$description

rownames(pp) <- traits

tab <- tab[, c("bp_beta", "sh_beta", "bp-sh_beta")]
tab <- as.matrix(tab)
rownames(tab) <- traits
colnames(tab) <- c("CBP PRS", "SGIT PRS", "CBP UGIT PRS")

#tab_p <- round(tab, digits = 2)
#tab_p[nonsig] <- NA

# Load cluster data
load("hclust_level_2_thr_1.25.RData")

# Select the most representative traits from clusters
table(names(h2) %in% rownames(tab))
which(!(rownames(tab) %in% names(h2)))
#150 151 152 153 154 155
rownames(tab)[c(150, 151, 152, 153, 154, 155)]
table(names(h2) == rownames(tab)[-c(150, 151, 152, 153, 154, 155)]) #all true
table(h2)
#h2
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
# 6  4 28 14 10  4  8  8 11  3  7  4 14  7  4  4  7  4  2
h2 <- sort(h2, decreasing = FALSE)
tab_cut <- tab[names(h2), ]
pp_cut <- pp[names(h2), ]

c1 <- which(h2 == 1)
min(pp_cut[c1, ])
table(pp_cut == min(pp_cut[c1, ])) # unique
i1 <- which(pp_cut == min(pp_cut[c1, ]))
traits <- c(rownames(which(pp_cut == pp_cut[i1], arr.ind = TRUE))) # trait from the 1st cluster

c2 <- which(h2 == 2)
min(pp_cut[c2, ])
table(pp_cut == min(pp_cut[c2, ])) # unique
i2 <- which(pp_cut == min(pp_cut[c2, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i2], arr.ind = TRUE))) # add trait from the 2nd cluster

c3 <- which(h2 == 3)
min(pp_cut[c3, ])
table(pp_cut == min(pp_cut[c3, ])) # unique
i3 <- which(pp_cut == min(pp_cut[c3, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i3], arr.ind = TRUE))) # add trait from the 3d cluster

c4 <- which(h2 == 4)
min(pp_cut[c4, ])
table(pp_cut == min(pp_cut[c4, ])) # unique
i4 <- which(pp_cut == min(pp_cut[c4, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i4], arr.ind = TRUE))) # add trait from the 4th cluster

c5 <- which(h2 == 5)
min(pp_cut[c5, ])
table(pp_cut == min(pp_cut[c5, ])) # unique
i5 <- which(pp_cut == min(pp_cut[c5, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i5], arr.ind = TRUE))) # add trait from the 5th cluster

c6 <- which(h2 == 6)
min(pp_cut[c6, ])
table(pp_cut == min(pp_cut[c6, ])) # unique
i6 <- which(pp_cut == min(pp_cut[c6, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i6], arr.ind = TRUE))) # add trait from the 6th cluster

c7 <- which(h2 == 7)
min(pp_cut[c7, ])
table(pp_cut == min(pp_cut[c7, ])) # unique
i7 <- which(pp_cut == min(pp_cut[c7, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i7], arr.ind = TRUE)))# add trait from the 7th cluster

c8 <- which(h2 == 8)
min(pp_cut[c8, ])
table(pp_cut == min(pp_cut[c8, ])) # unique
i8 <- which(pp_cut == min(pp_cut[c8, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i8], arr.ind = TRUE)))# add trait from the 8th cluster

c9 <- which(h2 == 9)
min(pp_cut[c9, ])
table(pp_cut == min(pp_cut[c9, ])) # unique
i9 <- which(pp_cut == min(pp_cut[c9, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i9], arr.ind = TRUE))) # add trait from the 9th cluster

c10 <- which(h2 == 10)
min(pp_cut[c10, ])
table(pp_cut == min(pp_cut[c10, ])) # unique
i10 <- which(pp_cut == min(pp_cut[c10, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i10], arr.ind = TRUE))) # add trait from the 10th cluster

c11 <- which(h2 == 11)
min(pp_cut[c11, ])
table(pp_cut == min(pp_cut[c11, ])) # unique
i11 <- which(pp_cut == min(pp_cut[c11, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i11], arr.ind = TRUE))) # add trait from the 11th cluster

c12 <- which(h2 == 12)
min(pp_cut[c12, ])
table(pp_cut == min(pp_cut[c12, ])) # unique
i12 <- which(pp_cut == min(pp_cut[c12, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i12], arr.ind = TRUE))) # add trait from the 12th cluster

c13 <- which(h2 == 13)
min(pp_cut[c13, ])
table(pp_cut == min(pp_cut[c13, ])) # unique
i13 <- which(pp_cut == min(pp_cut[c13, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i13], arr.ind = TRUE))) # add trait from the 13th cluster

c14 <- which(h2 == 14)
min(pp_cut[c14, ])
table(pp_cut == min(pp_cut[c14, ])) # unique
i14 <- which(pp_cut == min(pp_cut[c14, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i14], arr.ind = TRUE))) # add trait from the 14th cluster

c15 <- which(h2 == 15)
min(pp_cut[c15, ])
table(pp_cut == min(pp_cut[c15, ])) # unique
i15 <- which(pp_cut == min(pp_cut[c15, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i15], arr.ind = TRUE))) # add trait from the 15th cluster

c16 <- which(h2 == 16)
min(pp_cut[c16, ])
table(pp_cut == min(pp_cut[c16, ])) # unique
i16 <- which(pp_cut == min(pp_cut[c16, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i16], arr.ind = TRUE))) # add trait from the 16th cluster

c17 <- which(h2 == 17)
min(pp_cut[c17, ])
table(pp_cut == min(pp_cut[c17, ])) # unique
i17 <- which(pp_cut == min(pp_cut[c17, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i17], arr.ind = TRUE))) # add trait from the 17th cluster

c18 <- which(h2 == 18)
min(pp_cut[c18, ])
table(pp_cut == min(pp_cut[c18, ])) # unique
i18 <- which(pp_cut == min(pp_cut[c18, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i18], arr.ind = TRUE))) # add trait from the 18th cluster

c19 <- which(h2 == 19)
min(pp_cut[c19, ])
table(pp_cut == min(pp_cut[c19, ])) # unique
i19 <- which(pp_cut == min(pp_cut[c19, ]))
traits <- c(traits, rownames(which(pp_cut == pp_cut[i19], arr.ind = TRUE))) # add trait from the 19th cluster

h2[traits] # self-check, each cluster is represented
tab_cut <- tab_cut[traits, ]
pp_cut <- pp_cut[traits, ]

#clusters <- c("Pain-related traits and surgery", "Thrombosis", "Injuries", "Respiratory illnesses and smoking", "Anthropometric traits", "Osteoarthritis and other musculoskeletal disorders", "Socio-economic and family status", "Cardio-vascular traits", "Physical activity", "Psychometric traits", "Head pain and medication for pain relief")

#rownames(rg_m_cut) <- rownames(rg_p_cut) <- clusters

#thr <- 0.05 / (((324*324)-324)/2 + (730 - 322)*2)

tab_fin <- rbind(tab_cut, tab[c(150, 151, 152, 153, 154, 155), ])
pp_fin <- rbind(pp_cut, pp[c(150, 151, 152, 153, 154, 155), ])
tab_fin <- round(tab_fin, digits = 2)
notsig <- which(pp_fin > thr)
pp_fin[notsig] <- NA

pdf("heatmap_level_2_thr_1.25_test_nonrel.pdf", width = 50, height = 50)
heatmap.2(tab_fin, cellnote = tab_fin, notecex = 4.0, notecol = "black", na.color = "grey",
          col = redblue(75), margins = c(25, 120), density.info = "none", trace = "none",
          cexRow = 4.0, srtCol = 45, cexCol = 4.0,
          lhei = c(1, 8), lwid = c(1.2, 4), Colv = FALSE, dendrogram = "row")
dev.off()

save(traits, tab_fin, pp_fin, file = "heatmap_level_2_clusters_1.25.RData")
