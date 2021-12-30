# Aim of this script is to perform clustering of the ICD10 and OPS traits
# (Test sample, non-relatieves only; standardized PRS)

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/")

# Load tables of ICD10 and OPCS codes
load("./icd10_iid_cbp_prev_filtered_test_nonrelatives.RData")
load("./opcs_iid_cbp_prev_filtered_test_nonrelatives.RData")

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load glm() results
# BP PRS
load("./glm_icd_bp_prs_filtered_test_nonrelatives.RData")
load("./glm_opcs_bp_prs_filtered_test_nonrelatives.RData")

# SH PRS
load("./glm_icd_sh_prs_filtered_test_nonrelatives.RData")
load("./glm_opcs_sh_prs_filtered_test_nonrelatives.RData")

# BP-SH PRS
load("./glm_icd_bp_sh_prs_filtered_test_nonrelatives.RData")

ls()

# Select ICD10 and OPCS codes for clustering
codes <- names(bp_icd)
codes <- c(codes, names(bp_opcs))
codes <- c(codes, names(sh_icd))
codes <- c(codes, names(sh_opcs))
codes <- c(codes, names(bp_sh_icd))
length(codes)
length(unique(codes))
codes <- unique(codes)

# Calculate phenotypic correlation matrix between ICD10 and OPCS codes
icd_ind <- na.omit(match(codes, colnames(icd_f_prev_test_nonr)))
opcs_ind <- na.omit(match(codes, colnames(opcs_f_prev_test_nonr)))
pheno_joint <- cbind(icd_f_prev_test_nonr[, icd_ind], opcs_f_prev_test_nonr[, opcs_ind])
pheno_cor <- cor(pheno_joint)

# Perform hierarchical clustering
dd <- (1 - abs(pheno_cor)) # transform to distance matrix

d <- as.dist(dd)
h1 <- hclust(d, method = "ward.D2")
results <- data.table::fread("glm_results_extended_desc_test_nonrelatives.txt", data.table = F)
j <- match(h1$labels, results$code)
table(results$code[j] == h1$labels) # all true
h1$labels <- results$description[j]

pdf("hclust_test_nonrelatives.pdf", width = 30, height = 30)
plot(h1, hang = -1)
abline(h = 1.5, col = "blue")
abline(h = 1.4, col = "red")
abline(h = 1.3, col = "green")
dev.off()

h2 <- cutree(tree = h1, h = 1.5)
save(h2, file = "hclust_1.5.RData")

# Heatmap visualization
library(corrplot)

rg2plot <- pheno_cor[h1$order, h1$order]

pdf("heatmap_hclust_test_nonrelatieves.pdf", width = 30, height = 30)
corrplot(as.matrix(rg2plot), type = 'full', is.corr = T)
dev.off()




