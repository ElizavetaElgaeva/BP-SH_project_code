# Aim of this script is to create a table with glm() results for BP/SH/BP-SH PRS

setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Load full glm() results
# BP 
load("./icd10_vs_bp_prs.RData")
load("./opcs_vs_bp_prs.RData")

# SH
load("./icd10_vs_sh_prs.RData")
load("./opcs_vs_sh_prs.RData")

# BP-SH
load("./icd10_vs_bp_sh_prs.RData")
load("./opcs_vs_bp_sh_prs.RData")


# Load glm(0 results of interest
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

# Select ICD10 and OPCS codes of interest
codes <- names(bp_icd)
codes <- c(codes, names(bp_opcs))
codes <- c(codes, names(sh_icd))
codes <- c(codes, names(sh_opcs))
codes <- c(codes, names(icd_opcs))
length(codes)
length(unique(codes))
codes <- unique(codes)

icd_ind <- which(codes %in% names(bp_bin_cov))

# Create a table
tab <- as.data.frame(matrix(ncol = 15, nrow = length(codes)))
colnames(tab) <- c("code", "bp_or", "bp_beta", "bp_se", "bp_p", "sh_or", "sh_beta", "sh_se", "sh_p", "bp-sh_or", "bp-sh_beta", "bp-sh_se", "bp-sh_p", "nomenclature", "code_prevalence")
tab$code <- codes
tab$nomenclature[icd_ind] <- "ICD10"
tab$nomenclature[-icd_ind] <- "OPCS"

bp_full <- c(bp_bin_cov, opcs_bp_bin)
bp_full <- bp_full[codes]
tab$"bp_or" <- sapply(bp_full, function(x) exp(x[, 1]))
tab$"bp_beta" <- sapply(bp_full, function(x) x[, 1])
tab$"bp_se" <- sapply(bp_full, function(x) x[, 2])
tab$"bp_p" <- sapply(bp_full, function(x) x[, 4])

sh_full <- c(sh_bin_cov, opcs_sh_bin)
sh_full <- sh_full[codes]
tab$"sh_or" <- sapply(sh_full, function(x) exp(x[, 1]))
tab$"sh_beta" <- sapply(sh_full, function(x) x[, 1])
tab$"sh_se" <- sapply(sh_full, function(x) x[, 2])
tab$"sh_p" <- sapply(sh_full, function(x) x[, 4])

bp_sh_full <- c(bp_sh_bin_cov, opcs_bp_sh_bin)
bp_sh_full <- bp_sh_full[codes]
tab$"bp-sh_or" <- sapply(bp_sh_full, function(x) exp(x[, 1]))
tab$"bp-sh_beta" <- sapply(bp_sh_full, function(x) x[, 1])
tab$"bp-sh_se" <- sapply(bp_sh_full, function(x) x[, 2])
tab$"bp-sh_p" <- sapply(bp_sh_full, function(x) x[, 4])

load("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/prev_of_icd10_filtered.RData")
load("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/prev_opcs_filtered.RData")
prev_full <- c(prev_icd_f, prev_opcs_f)
prev_full <- prev_full[codes]
tab$"code_prevalence" <- prev_full

data.table::fwrite(tab, file = "glm_results.txt", sep = "\t", dec = ".")

# Extend the table
load("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/bp_prs_bin_iid_filtered.RData")
pain <- data.table::fread("/home/freydin/Pheno/mv_mar18_chron.phen.txt", data.table = F)
ind <- match(bp_f_bin$IID, pain$IID)
pain <- pain[ind, ]
pain <- pain[, c("back", "neck", "hip", "stom", "knee", "head")]
pain_bp_cov <- lapply(pain, function(x) 
       		    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp_bin, data = bp_f_bin, family = "binomial") })
pain_bp_cov <- lapply(pain_bp_cov, function(x) tail(summary(x)$coefficients, n = 1))

pain_sh_cov <- lapply(pain, function(x)
                    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_sh_bin, data = bp_f_bin, family = "binomial") })
pain_sh_cov <- lapply(pain_sh_cov, function(x) tail(summary(x)$coefficients, n = 1))

pain_bp_sh_cov <- lapply(pain, function(x)
                    { glm(x ~ Age + Sex + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + prs_bp_sh_bin, data = bp_f_bin, family = "binomial") })
pain_bp_sh_cov <- lapply(pain_bp_sh_cov, function(x) tail(summary(x)$coefficients, n = 1))

tab2 <- as.data.frame(matrix(ncol = 15, nrow = ncol(pain)))
colnames(tab2) <- c("code", "bp_or", "bp_beta", "bp_se", "bp_p", "sh_or", "sh_beta", "sh_se", "sh_p", "bp-sh_or", "bp-sh_beta", "bp-sh_se", "bp-sh_p", "nomenclature", "code_prevalence")
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

data.table::fwrite(tab_joint, file = "glm_results_extended.txt", sep = "\t", dec = ".")

# Add description
icd_des <- readxl::read_excel("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/ICD10_coding19.xlsx", sheet = 1)
opcs_des <- readxl::read_excel("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/OPCS_coding240.xlsx", sheet = 1)

i_icd_des <- match(codes[icd_ind], icd_des$"coding")
i_opcs_des <- match(codes[-icd_ind], opcs_des$"coding")

icd_des <- icd_des[i_icd_des, c("coding", "meaning")]
opcs_des <- opcs_des[i_opcs_des, c("coding", "meaning")]
icd_opcs_des <- rbind(icd_des, opcs_des)
i2 <- match(tab_joint$code, icd_opcs_des$"coding") 

tab_joint$description <- icd_opcs_des$meaning[i2]

data.table::fwrite(tab_joint, file = "glm_results_extended_desc.txt", sep = "\t", dec = ".")
