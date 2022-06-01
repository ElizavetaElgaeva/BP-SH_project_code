# Aim of this script is to create a table with glm() results for BP/SH/BP-SH PRS
# (Test sample, non-relatieves only; standardized PRS; ICD10 and OPCS codes combined at level 2)

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load full glm() results
# BP 
load("./icd10_level_2_vs_bp_prs_test_nonrelatives_var1.RData")
load("./opcs_level_2_vs_bp_prs_test_nonrelatives_var1.RData")

# SH
load("./icd10_level_2_vs_sh_prs_test_nonrelatives_var1.RData")
load("./opcs_level_2_vs_sh_prs_test_nonrelatives_var1.RData")

# BP-SH
load("./icd10_level_2_vs_bp_sh_prs_test_nonrelatives_var1.RData")
load("./opcs_level_2_vs_bp_sh_prs_test_nonrelatives_var1.RData")


# Load glm results filtered by statistical significance
# BP PRS
load("./glm_icd_level_2_bp_prs_filtered_test_nonrelatives.RData")
load("./glm_opcs_level_2_bp_prs_filtered_test_nonrelatives.RData")

# SH PRS
load("./glm_icd_level_2_sh_prs_filtered_test_nonrelatives.RData")
load("./glm_opcs_level_2_sh_prs_filtered_test_nonrelatives.RData")

# BP-SH PRS
load("./glm_opcs_level_2_bp_sh_prs_filtered_test_nonrelatives.RData")

ls()

# Select ICD10 and OPCS codes of interest
codes <- names(bp_icd)
codes <- c(codes, names(sh_icd))
codes <- unique(codes)

codes2 <- names(bp_opcs)
codes2 <- c(codes2, names(sh_opcs))
codes2 <- c(codes2, names(bp_sh_opcs))
codes2 <- unique(codes2)

intersect(codes, codes2) # "E03" "K63" "M47" "Z82" "Z85" "Z86" "Z90" "Z92" "J18"
all_codes <- c(codes, codes2) # length(codes) = 135, length(codes2) = 87
icd_ind <- c(1:length(codes))

# Create a table
tab <- as.data.frame(matrix(ncol = 15, nrow = length(all_codes)))
colnames(tab) <- c("nomenclature", "code", "bp_or", "bp_beta", "bp_se", "bp_p", "sh_or", "sh_beta", "sh_se", "sh_p", "bp-sh_or", "bp-sh_beta", "bp-sh_se", "bp-sh_p", "code_prevalence")
tab$code <- all_codes
tab$nomenclature[icd_ind] <- "ICD10"
tab$nomenclature[-icd_ind] <- "OPCS"

bp_full <- c(bp_test_nonr_var1_cov[codes], opcs_bp_test_nonr_var1[codes2])
tab$"bp_or" <- sapply(bp_full, function(x) exp(x[, 1]))
tab$"bp_beta" <- sapply(bp_full, function(x) x[, 1])
tab$"bp_se" <- sapply(bp_full, function(x) x[, 2])
tab$"bp_p" <- sapply(bp_full, function(x) x[, 4])

sh_full <- c(sh_test_nonr_var1_cov[codes], opcs_sh_test_nonr_var1[codes2])
tab$"sh_or" <- sapply(sh_full, function(x) exp(x[, 1]))
tab$"sh_beta" <- sapply(sh_full, function(x) x[, 1])
tab$"sh_se" <- sapply(sh_full, function(x) x[, 2])
tab$"sh_p" <- sapply(sh_full, function(x) x[, 4])

bp_sh_full <- c(bp_sh_test_nonr_var1_cov[codes], opcs_bp_sh_test_nonr_var1[codes2])
tab$"bp-sh_or" <- sapply(bp_sh_full, function(x) exp(x[, 1]))
tab$"bp-sh_beta" <- sapply(bp_sh_full, function(x) x[, 1])
tab$"bp-sh_se" <- sapply(bp_sh_full, function(x) x[, 2])
tab$"bp-sh_p" <- sapply(bp_sh_full, function(x) x[, 4])

load("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/prev_of_icd10_level_2_filtered.RData")
load("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/prev_opcs_level_2_filtered.RData")
prev_full <- c(prev_icd_f[codes], prev_opcs_f[codes2])
tab$"code_prevalence" <- prev_full

data.table::fwrite(tab, file = "glm_results_level_2_codes_test_nonrel.txt", sep = "\t", dec = ".")

# Extend the table
load("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/bp_prs_iid_icd10_filtered_test_nonrelatives_prs_var1.RData")
pain <- data.table::fread("/home/freydin/Pheno/mv_mar18_chron.phen.txt", data.table = F)
ind <- match(bp_f_test_nonr_var1$IID, pain$IID)
pain <- pain[ind, ]
pain <- pain[, c("back", "neck", "hip", "stom", "knee", "head")]
pain_bp_cov <- lapply(pain, function(x) 
       		    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp, data = bp_f_test_nonr_var1, family = "binomial") })
pain_bp_cov <- lapply(pain_bp_cov, function(x) tail(summary(x)$coefficients, n = 1))

pain_sh_cov <- lapply(pain, function(x)
                    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_sh, data = bp_f_test_nonr_var1, family = "binomial") })
pain_sh_cov <- lapply(pain_sh_cov, function(x) tail(summary(x)$coefficients, n = 1))

pain_bp_sh_cov <- lapply(pain, function(x)
                    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp_sh, data = bp_f_test_nonr_var1, family = "binomial") })
pain_bp_sh_cov <- lapply(pain_bp_sh_cov, function(x) tail(summary(x)$coefficients, n = 1))

tab2 <- as.data.frame(matrix(ncol = 15, nrow = ncol(pain)))
colnames(tab2) <- c("nomenclature", "code", "bp_or", "bp_beta", "bp_se", "bp_p", "sh_or", "sh_beta", "sh_se", "sh_p", "bp-sh_or", "bp-sh_beta", "bp-sh_se", "bp-sh_p", "code_prevalence")
tab2$code <- colnames(pain)

tab2$"bp_or" <- sapply(pain_bp_cov, function(x) exp(x[, 1]))
tab2$"bp_beta" <- sapply(pain_bp_cov, function(x) x[, 1])
tab2$"bp_se" <- sapply(pain_bp_cov, function(x) x[, 2])
tab2$"bp_p" <- sapply(pain_bp_cov, function(x) x[, 4])

tab2$"sh_or" <- sapply(pain_sh_cov, function(x) exp(x[, 1]))
tab2$"sh_beta" <- sapply(pain_sh_cov, function(x) x[, 1])
tab2$"sh_se" <- sapply(pain_sh_cov, function(x) x[, 2])
tab2$"sh_p" <- sapply(pain_sh_cov, function(x) x[, 4])

tab2$"bp-sh_or" <- sapply(pain_bp_sh_cov, function(x) exp(x[, 1]))
tab2$"bp-sh_beta" <- sapply(pain_bp_sh_cov, function(x) x[, 1])
tab2$"bp-sh_se" <- sapply(pain_bp_sh_cov, function(x) x[, 2])
tab2$"bp-sh_p" <- sapply(pain_bp_sh_cov, function(x) x[, 4])

n <- nrow(pain)
pain_prev <-  lapply(pain, function(x) sum(x)/n * 100)
tab2$"code_prevalence" <- unlist(pain_prev)

tab_joint <- rbind(tab, tab2)

data.table::fwrite(tab_joint, file = "glm_results_level_2_codes_extended_test_nonrelatives.txt", sep = "\t", dec = ".")

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
data.table::fwrite(tab_joint, file = "glm_results_level_2_extended_desc_test_nonrelatives.txt", sep = "\t", dec = ".")
