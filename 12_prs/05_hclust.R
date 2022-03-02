# Aim of this script is to perform clustering of the ICD10 and OPCS traits

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/")

# Load tables of ICD10 and OPCS codes
load("./icd10_iid_cbp_prev_filtered.RData")
load("./opcs_iid_cbp_prev_filtered.RData")

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load glm() results
# BP PRS
load("./glm_icd_bp_prs_filtered.RData")
load("./glm_opcs_bp_prs_filtered.RData")

# SH PRS
load("./glm_icd_sh_prs_filtered.RData")
load("./glm_opcs_sh_prs_filtered.RData")

# BP-SH PRS
load("./glm_icd_bp_sh_prs_filtered.RData")
load("./glm_opcs_bp_sh_prs_filtered.RData")
load("./glm_icd_opcs_bp_sh_prs_filtered.RData")

ls()

# Select ICD10 and OPCS codes for clustering
codes <- names(bp_icd)
codes <- c(codes, names(bp_opcs))
codes <- c(codes, names(sh_icd))
codes <- c(codes, names(sh_opcs))
codes <- c(codes, names(icd_opcs))
length(codes)
length(unique(codes))
codes <- unique(codes)

# Calculate phenotypic correlation matrix between ICD10 and OPCS codes
icd_ind <- na.omit(match(codes, colnames(icd_f_prev)))
opcs_ind <- na.omit(match(codes, colnames(opcs_f_prev)))
pheno_joint <- cbind(icd_f_prev[, icd_ind], opcs_f_prev[, opcs_ind])
pheno_cor <- cor(pheno_joint)

# Perform hierarchical clustering
dd <- (1 - abs(pheno_cor)) # transform to distance matrix

d <- as.dist(dd)
h1 <- hclust(d, method = "ward.D2")

pdf("hclust.pdf", width = 30, height = 30)
plot(h1, hang = -1)
dev.off()

# Heatmap visualization
library(corrplot)

rg2plot <- pheno_cor[h1$order, h1$order]

pdf("heatmap_hclust.pdf", width = 30, height = 30)
corrplot(as.matrix(rg2plot), type = 'full', is.corr = T)
dev.off()




