# Aim of this script is to select non-relatieves from the test sampe for PRSanalysis
# and filter codes by prevalence

library(data.table)
setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Read file with non-relatives IIDs
nonrel <- fread("ukbb.nonrel.txt", data.table = F)
dim(nonrel)
head(nonrel)

# Read file with test sample IIDs
test <- fread("test_id.txt", data.table = F)
dim(test)
head(test)

# Find the intersect between list of non-relatieves and test sample
inter <- intersect(nonrel$V1, test$V1)
length(inter) # 120217

# Let's extract the IIDs of interest from tables with PRS, ICD10 and OPCS data
setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

load("./bp_prs_iid_icd10_filtered.RData") # BP data
load("./icd10_iid_level_2_cbp_filtered.RData") # ICD10 data
load("./opcs_level_2_iid_cbp_filtered.RData") # OPCS data

ls()

table(inter %in% bp_f$IID)
table(inter %in% icd_f_l2$IID)
table(inter %in% opcs_f_l2$IID)

# IID order is the same in bp_f, icd_f_l2, opcs_f_l2
i <- which(bp_f$IID %in% inter)
bp_f_test_nonr <- bp_f[i, ]
icd_f_l2_test_nonr <- icd_f_l2[i, ]
opcs_f_l2_test_nonr <- opcs_f_l2[i, ]

# Filter ICD10 codes by prevalence
prev <- lapply(icd_f_l2_test_nonr[, -c(1:22)], function(x) sum(x)) # the first 22 columns contain no ICD10 codes
prev <- unlist(prev)
n <- nrow(icd_f_l2_test_nonr)
prev <- lapply(prev, function(x) x/n*100)
ind_prev <- which(prev > 0.5 & prev < 99.5) # 245 ICD10 codes meet the criteria
icd_f_l2_test_nonr_prev <- icd_f_l2_test_nonr[, c(1:22, ind_prev + 22)] # need to shift the indexes since the were defined in data without first 22 columns
prev_icd_f <- prev[ind_prev]

# Filter OPCS codes by prevalence
prev <- lapply(opcs_f_l2_test_nonr[, -c(1:22)], function(x) sum(x)) # the first 22 columns contain no OPCS codes
prev <- unlist(prev)
n <- nrow(opcs_f_l2_test_nonr)
prev <- lapply(prev, function(x) x/n*100)
ind_prev <- which(prev > 0.5 & prev < 99.5) # 201 OPCS codes meet the criteria
opcs_f_l2_test_nonr_prev <- opcs_f_l2_test_nonr[, c(1:22, ind_prev + 22)] # need to shift the indexes since the were defined in data without first 22 columns
prev_opcs_f <- prev[ind_prev]

# Save data
save(bp_f_test_nonr, file = "bp_prs_iid_icd10_filtered_test_nonrelatives.RData")
save(icd_f_l2_test_nonr_prev, file = "icd10_iid_level_2_cbp_prev_filtered_test_nonrelatives.RData")
save(opcs_f_l2_test_nonr_prev, file = "opcs_iid_level_2_cbp_prev_filtered_test_nonrelatives.RData")
save(prev_icd_f, file = "prev_of_icd10_level_2_filtered.RData")
save(prev_opcs_f, file = "prev_opcs_level_2_filtered.RData")
