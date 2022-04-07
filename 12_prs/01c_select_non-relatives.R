# Aim of this script is to select non-relatieves from the test sampe for PRSanalysis

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
load("./icd10_iid_cbp_prev_filtered.RData") # ICD10 data
load("./opcs_iid_cbp_prev_filtered.RData") # OPCS data

ls()

table(inter %in% bp_f$IID)
table(inter %in% icd_f_prev$IID)
table(inter %in% opcs_f_prev$IID)

# IID order is the same in bp_f, icd_f_prev, opcs_f_prev
i <- which(bp_f$IID %in% inter)
bp_f_test_nonr <- bp_f[i, ]
icd_f_prev_test_nonr <- icd_f_prev[i, ]
opcs_f_prev_test_nonr <- opcs_f_prev[i, ]

save(bp_f_test_nonr, file = "bp_prs_iid_icd10_filtered_test_nonrelatives.RData")
save(icd_f_prev_test_nonr, file = "icd10_iid_cbp_prev_filtered_test_nonrelatives.RData")
save(opcs_f_prev_test_nonr, file = "opcs_iid_cbp_prev_filtered_test_nonrelatives.RData")

