# Aim of this script is to filter OPCS data

library(data.table)

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/")

# Read table with ukbb data OPCS codes 
opcs <- fread("./OPCS_main_secondary_merged.txt", data.table = F)

# Read clean CBP PRS tableload("./bp_prs_iid_icd10_filtered.RData")
load("./bp_prs_iid_icd10_filtered.RData")

# Filter OPCS data to keep only those IIDs, that are in clean BP table
nrow(opcs)
#[1] 487322
nrow(bp_f)
#[1] 439762
table(opcs$IID %in% bp_f$IID)
# FALSE   TRUE
# 47560 439762
table(bp_f$IID %in% opcs$IID)
# TRUE
# 439762

ind <- match(bp_f$IID, opcs$IID)
opcs_f <- opcs[ind, ]

table(opcs_f$IID %in% bp_f$IID)

# Check if there non-binary traits
# columns 1 and the last 21 are not OPCS codes
bin <- lapply(opcs_f, function(x) { all(x %in% 0:1) })
table(unlist(bin)[-c(1, 8003:8023)] == TRUE)
# TRUE
#8001
# no non-binary traits

# Check for NAs
na <- lapply(opcs_f, function(x) all(!is.na(x)))
table(unlist(na)[-c(8003:8023)] == TRUE)
# TRUE
#8002

dim(na.omit(opcs_f))
#[1] 437604  8023

# Check if empty
empty <- lapply(opcs_f, function(x) all(x != ""))
table(unlist(empty)[-c(8003:8023)] == TRUE)
# TRUE
#8002
# no emty fields

# Filter by prevalence
prev <- lapply(opcs_f[, -c(1, 8003:8023)], function(x) sum(x))
prev <- unlist(prev)
n <- nrow(opcs_f)
prev <- lapply(prev, function(x) x/n*100)
ind_prev <- which(prev > 0.5 & prev < 99.5)
opcs_f_prev <- opcs_f[, c(1, ind_prev + 1)]

# Save data
opcs_f <- opcs_f[, -c(8003:8023)]
readme_opcs_f <- "OPCS table with only IIDs presented in bp_prs_iid_icd10_filtered.RData"
save(opcs_f, readme_opcs_f, file = "./opcs_iid_cbp_filtered.RData")

readme_opcs_f_prev <- "OPCS table with only IIDs presented in bp_prs_iid_icd10_filtered.RData and OPCS codes for binary traits with prevalence > 0.5% and < 99.5%"
save(opcs_f_prev, readme_opcs_f_prev, file = "./opcs_iid_cbp_prev_filtered.RData")


