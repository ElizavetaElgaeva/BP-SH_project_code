# Aim of this script is to create a table with glm() results for BP/SH/BP-SH PRS
# (Test sample, non-relatieves only; standardized PRS; ICD10 and OPCS codes combined at level 2)
# Note: ICD10 from chapter I - XVII only; X, Y, Z OPCS codes excluded

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load full glm() results
# BP 
load("./icd10_level_2_chapter_1-17_vs_bp_prs_test_nonrelatives_var1.RData")
load("./opcs_level_2_no_xyz_vs_bp_prs_test_nonrelatives_var1.RData")

# Load glm results filtered by statistical significance
# BP PRS
load("./glm_icd_level_2_chapter_1-17_bp_prs_filtered_test_nonrelatives_for_PRS_paper.RData")
load("./glm_opcs_level_2_no_xyz_bp_prs_filtered_test_nonrelatives_for_PRS_paper.RData")

ls()

# Select ICD10 and OPCS codes
codes <- names(bp_test_nonr_var1_cov)

codes2 <- names(opcs_bp_test_nonr_var1)

intersect(codes, codes2) # "C79" "E03" "E14" "F10" "G45" "H25" "H33" "J18" "K40" "K60" "K63" "M47" "M65" "N17" "N30"
all_codes <- c(codes, codes2) # length(codes) = 165, length(codes2) = 132
icd_ind <- c(1:length(codes))

# Create a table
tab <- as.data.frame(matrix(ncol = 7, nrow = length(all_codes)))
colnames(tab) <- c("nomenclature", "code", "bp_or", "bp_beta", "bp_se", "bp_p", "code_prevalence")
tab$code <- all_codes
tab$nomenclature[icd_ind] <- "ICD10"
tab$nomenclature[-icd_ind] <- "OPCS"

bp_full <- c(bp_test_nonr_var1_cov, opcs_bp_test_nonr_var1)
tab$"bp_or" <- sapply(bp_full, function(x) exp(x[, 1]))
tab$"bp_beta" <- sapply(bp_full, function(x) x[, 1])
tab$"bp_se" <- sapply(bp_full, function(x) x[, 2])
tab$"bp_p" <- sapply(bp_full, function(x) x[, 4])

load("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/prev_of_icd10_level_2_filtered.RData")
load("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/prev_opcs_level_2_filtered.RData")
prev_full <- c(prev_icd_f[codes], prev_opcs_f[codes2])
tab$"code_prevalence" <- prev_full

data.table::fwrite(tab, file = "glm_results_level_2_preselected_codes_test_nonrel_for_PRS_paper.txt", sep = "\t", dec = ".")

# Extend the table
load("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/bp_prs_iid_icd10_filtered_test_nonrelatives_prs_var1.RData")
pain <- data.table::fread("/home/freydin/Pheno/mv_mar18_chron.phen.txt", data.table = F)
ind <- match(bp_f_test_nonr_var1$IID, pain$IID)
pain <- pain[ind, ]
pain <- pain[, c("back", "neck", "hip", "stom", "knee", "head")]
pain_bp_cov <- lapply(pain, function(x) 
       		    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp, data = bp_f_test_nonr_var1, family = "binomial") })
pain_bp_cov <- lapply(pain_bp_cov, function(x) tail(summary(x)$coefficients, n = 1))

tab2 <- as.data.frame(matrix(ncol = 7, nrow = ncol(pain)))
colnames(tab2) <- c("nomenclature", "code", "bp_or", "bp_beta", "bp_se", "bp_p", "code_prevalence")
tab2$code <- colnames(pain)

tab2$"bp_or" <- sapply(pain_bp_cov, function(x) exp(x[, 1]))
tab2$"bp_beta" <- sapply(pain_bp_cov, function(x) x[, 1])
tab2$"bp_se" <- sapply(pain_bp_cov, function(x) x[, 2])
tab2$"bp_p" <- sapply(pain_bp_cov, function(x) x[, 4])

n <- nrow(pain)
pain_prev <-  lapply(pain, function(x) sum(x)/n * 100)
tab2$"code_prevalence" <- unlist(pain_prev)

tab_joint <- rbind(tab, tab2)

data.table::fwrite(tab_joint, file = "glm_results_level_2_preselected_codes_extended_test_nonrelatives_for_PRS_paper.txt", sep = "\t", dec = ".")

# Add description
icd_des <- readxl::read_excel("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/ICD10_coding19.xlsx", sheet = 1)
opcs_des <- readxl::read_excel("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/OPCS_coding240.xlsx", sheet = 1)

i_icd_des <- match(codes, icd_des$"coding")
i_opcs_des <- match(codes2, opcs_des$"coding")

icd_des <- icd_des[i_icd_des, c("coding", "meaning")]
opcs_des <- opcs_des[i_opcs_des, c("coding", "meaning")]
icd_opcs_des <- rbind(icd_des, opcs_des)
i2 <- c(1:length(all_codes)) 
table(tab_joint$code[i2] == icd_opcs_des$"coding")
tab_joint$description <- NA
tab_joint$description[i2] <- icd_opcs_des$meaning
tab_joint$description[-i2] <- c("Chronic back pain", "Chronic neck pain", "Chronic hip pain", "Chronic stomach pain", "Chronic knee pain", "Chronic headache")
data.table::fwrite(tab_joint, file = "glm_results_level_2_preselected_codes_extended_desc_test_nonrelatives_for_PRS_paper.txt", sep = "\t", dec = ".")
