# Aim of this script is to create a heatmap visualizing glm() results

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

library(gplots)
library(data.table)

tab <- fread("glm_results_extended_desc.txt")
colnames(tab)

pp <- tab[, c("bp_p", "sh_p", "bp-sh_p")]
thr <- 3.646973e-05
nonsig <- which(pp > thr)

traits <- tab$"description"
traits[20:25] <- c("Chronic back pain", "Chronic neck pain", "Chronic hip pain", "Chronic stomach pain", "Chronic knee pain", "Chronic headache")
tab <- tab[, c("bp_beta", "sh_beta", "bp-sh_beta")]
tab <- as.matrix(tab)
rownames(tab) <- traits
colnames(tab) <- c("CBP PRS", "SGIT PRS", "CBP UGIT PRS")

tab_p <- round(tab, digits = 2)
tab_p[nonsig] <- NA

pdf("heatmap.pdf", width = 50, height = 50)
heatmap.2(tab, cellnote = tab_p, notecex = 4.0, notecol = "black",
	  col = redblue(75), margins = c(25, 120), density.info = "none", trace = "none",
	  cexRow = 4.0, srtCol = 45, cexCol = 4.0,
	  lhei = c(1, 8), lwid = c(1.2, 4), Colv = FALSE, dendrogram = "row")
dev.off()
