# Aim of this script is to visualize matrices of genetic annd phenotypic correlations
# for original traits, sh and traits minus sh

library(corrplot)
library(data.table)

setwd('/mnt/polyomica/projects/bp-sh/data/09_heatmaps')

# Specify p-value threshold
thr <- 0.05/28

# Load genetic correlations matrix
gcor <- read.table('./gene_corr_matrix.txt')
colnames(gcor) <- rownames(gcor)

h2_table <- read.csv('./gene_corr/h2.csv')

# Form a p-value matrix
cor_files <- list.files('./gene_corr', full.names = T, pattern = 'gene_corr_.+csv')
cor <- lapply(cor_files, fread)
dim(cor[[1]])
colnames(cor[[1]])

p_matrix <- as.data.table(matrix(NA, nrow = length(cor_files), ncol = length(cor_files)))
colnames(p_matrix) <- rownames(p_matrix) <- rownames(gcor)

for(i in 1:length(cor_files)){

	index <- match(rownames(gcor), cor[[i]]$"gwas_id_2")
	cor[[i]] <- cor[[i]][index, ]
	p_matrix[, i] <- cor[[i]]$pval

}
p_matrix <- as.matrix(p_matrix)
p_matrix[upper.tri(p_matrix, diag = F)] <- NA
p_matrix <- Matrix::forceSymmetric(p_matrix, uplo = "L")
p_matrix <- as.matrix(p_matrix)
p_matrix <- as.data.table(p_matrix)
colnames(p_matrix) <- rownames(p_matrix) <- rownames(gcor)

fwrite(p_matrix, 'gene_cor_p_val_matrix.txt', row.names = T, col.names = T, quote = F, sep = "\t", dec = ".")

p_matrix <- as.data.frame(p_matrix)

# Heretability vector reordering
ind_h2 <- match(rownames(gcor), h2_table$gwas_id)
h2 <- h2_table$h2[ind_h2] # obtain heritabilities

# Rename 
traits <- c("Hip", "Back", "Neck", "Knee", "Head", "Stomach", "SGC", "UGC of CBP")
colnames(gcor) <- rownames(gcor) <- traits
colnames(p_matrix) <- rownames(p_matrix) <- traits

# Heatmap rg plot
gcor <- as.matrix(gcor)
diag(gcor) <- h2
p_matrix <- as.matrix(p_matrix)

out <- 'heatmap_rg.pdf'
pdf(out, height = 7, width = 7)
corrplot(gcor, method = "square", tl.col = "black", p.mat = p_matrix, sig.level = thr, addCoef.col = "black", cl.cex=1)
dev.off()

# Load phenotypic correlations matrix
pheno <- read.table("./pheno_corr_matrix.txt")
colnames(pheno) <- rownames(pheno) <- traits
pheno <- as.matrix(pheno)
pheno[upper.tri(pheno, diag = F)] <- NA
pheno <- Matrix::forceSymmetric(pheno, uplo = "L")
pheno <- as.matrix(pheno)


# Heatmap phenotypic correlations plot
out <- 'heatmap_pheno.pdf'
pdf(out, height = 7, width = 7)
corrplot(pheno, is.corr = T, method = "square", tl.col = "black", addCoef.col = "black", cl.cex=1, cl.lim = c(-1, 1))
dev.off()



