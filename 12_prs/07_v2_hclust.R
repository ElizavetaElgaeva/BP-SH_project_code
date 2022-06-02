# Aim of this script is to perform clustering of the ICD10 and OPS traits
# (Test sample, non-relatieves only; standardized PRS; codes combined to level 2;
# ICD10 codes from chapters I-XVII only, OPCS codes without X, Y, Z chapter)

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/")

# Load tables of ICD10 and OPCS codes
load("./icd10_iid_level_2_cbp_prev_filtered_test_nonrelatives.RData")
load("./opcs_iid_level_2_cbp_prev_filtered_test_nonrelatives.RData")

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load glm() results
# BP PRS
load("./glm_icd_level_2_chapter_1-17_bp_prs_filtered_test_nonrelatives.RData")
load("./glm_opcs_level_2_no_xyz_bp_prs_filtered_test_nonrelatives.RData")

# SH PRS
load("./glm_icd_level_2_chapter_1-17_sh_prs_filtered_test_nonrelatives.RData")
load("./glm_opcs_level_2_no_xyz_sh_prs_filtered_test_nonrelatives.RData")

# BP-SH PRS
load("./glm_opcs_level_2_no_xyz_bp_sh_prs_filtered_test_nonrelatives.RData")

ls()

# Select ICD10 and OPCS codes for clustering
codes <- names(bp_icd)
codes <- c(codes, names(sh_icd))
codes <- unique(codes)

codes2 <- names(bp_opcs)
codes2 <- c(codes2, names(sh_opcs))
codes2 <- c(codes2, names(bp_sh_opcs))
codes2 <- unique(codes2)

intersect(codes, codes2) # "E03" "K63" "M47" "J18"
all_codes <- c(codes, codes2)

# Calculate phenotypic correlation matrix between ICD10 and OPCS codes
icd_ind <- na.omit(match(codes, colnames(icd_f_l2_test_nonr_prev)))
opcs_ind <- na.omit(match(codes2, colnames(opcs_f_l2_test_nonr_prev)))
pheno_joint <- cbind(icd_f_l2_test_nonr_prev[, icd_ind], opcs_f_l2_test_nonr_prev[, opcs_ind])
pheno_cor <- cor(pheno_joint)

# Perform hierarchical clustering
dd <- (1 - abs(pheno_cor)) # transform to distance matrix

colnames(dd)[1:92] <- rownames(dd)[1:92] <- paste0(colnames(dd)[1:92], "_ICD10")
colnames(dd)[93:149] <- rownames(dd)[93:149] <- paste0(colnames(dd)[93:149], "_OPCS4")

d <- as.dist(dd)
h1 <- hclust(d, method = "ward.D2")
results <- data.table::fread("glm_results_level_2_preselected_codes_extended_desc_test_nonrelatives.txt", data.table = F)
table(substr(h1$labels, 1, 3) %in% results$code) # all true
table(substr(h1$labels, 1, 3) == results$code[1:149]) # all true
j <- c(1:length(h1$labels))
table(results$code[j] == substr(h1$labels, 1, 3)) # all true
h1$labels <- results$description[j]

pdf("hclust_level_2_test_nonrelatives.pdf", width = 30, height = 30)
plot(h1, hang = -1)
abline(h = 1.5, col = "blue")
abline(h = 1.35, col = "red")
abline(h = 1.25, col = "green")
dev.off()

h2 <- cutree(tree = h1, h = 1.25)
save(h2, file = "hclust_level_2_thr_1.25.RData")

h3 <- cutree(tree = h1, h = 1.35)
save(h3, file = "hclust_level_2_thr_1.35.RData")

# Heatmap visualization
library(corrplot)

colnames(pheno_cor) <- rownames(pheno_cor) <- colnames(dd)
rg2plot <- pheno_cor[h1$order, h1$order]

pdf("heatmap_hclust_lenel_2_test_nonrelatieves.pdf", width = 30, height = 30)
corrplot(as.matrix(rg2plot), type = 'full', is.corr = T)
dev.off()




