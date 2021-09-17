# Aim of this script is to create a heatmap visualizing glm() results

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

library(gplots)
library(data.table)

tab <- fread("glm_results_extended_desc.txt")
colnames(tab)

pp <- tab[, c("bp_p", "sh_p", "bp-sh_p")]

traits <- tab$"code"
traits[20:25] <- c("Chronic back pain", "Chronic neck pain", "Chronic hip pain", "Chronic stomach pain", "Chronic knee pain", "Chronic headache")
tab <- tab[, c("bp_beta", "sh_beta", "bp-sh_beta")]
tab <- as.matrix(tab)
rownames(tab) <- traits
colnames(tab) <- c("CBP PRS", "SGIT PRS", "CBP UGIT PRS")

pdf("heatmap.pdf", width = 30, height = 30)
heatmap.2(as.matrix(tab), cellnote = round(tab, digits = 2),
	  notecol = "black", col = redblue(75),
	  margins = c(25,5), density.info = "none",
	  trace = "none", cexRow = 0.7,
	  srtCol = 45, cexCol = 0.6)
dev.off()
