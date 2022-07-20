# Aim of this script is to select non-relatieves; cases only; from the whole sampe of Europeans
# and filter codes by prevalence

library(data.table)
setwd("/mnt/polyomica/projects/bp-sh/data/12_prs")

# Read file with non-relatives IIDs
nonrel <- fread("ukbb.nonrel.txt", data.table = F)
dim(nonrel)
head(nonrel)

# Let's extract the IIDs of interest from tables with PRS, ICD10 and OPCS data
setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

load("./bp_prs_iid_icd10_filtered.RData") # BP data
load("icd10_iid_cbp_filtered.RData") # ICD10 data; not-combined codes
load("opcs_iid_cbp_filtered.RData") # OPCS data; not-combined codes

ls()

# Order IIDs

table(nonrel$V1 %in% bp_f$IID)
i1 <- which(nonrel$V1 %in% bp_f$IID)
nonrel <- nonrel[i1, ]
i2 <- match(nonrel, bp_f$IID)
bp_f <- bp_f[i2, ]
table(nonrel == bp_f$IID)

table(nonrel == icd_f$IID[i2])
icd_f <- icd_f[i2, ]

table(nonrel == opcs_f$IID[i2])
opcs_f <- opcs_f[i2, ]

# Select back pain cases only
table(bp_f$back)
#     0      1
#299970  65011
cases <- which(bp_f$back == 1)
bp_cases_nonr <- bp_f[cases, ]
icd_cases_nonr <- icd_f[cases, ]
opcs_cases_nonr <- opcs_f[cases, ]


# Filter ICD10 codes by prevalence
prev <- lapply(icd_cases_nonr[, -c(1)], function(x) sum(x)) 
prev <- unlist(prev)
n <- nrow(icd_cases_nonr)
prev <- lapply(prev, function(x) x/n*100)
ind_prev <- which(prev > 0.5 & prev < 99.5) # 375 ICD10 codes meet the criteria
icd_cases_nonr_prev <- icd_cases_nonr[, c(1, ind_prev + 1)]
prev_icd_cases_f <- prev[ind_prev]

# Filter OPCS codes by prevalence
prev <- lapply(opcs_cases_nonr[, -c(1)], function(x) sum(x))
prev <- unlist(prev)
n <- nrow(opcs_cases_nonr)
prev <- lapply(prev, function(x) x/n*100)
ind_prev <- which(prev > 0.5 & prev < 99.5) # 265 OPCS codes meet the criteria
opcs_cases_nonr_prev <- opcs_cases_nonr[, c(1, ind_prev + 1)]
prev_opcs_cases_f <- prev[ind_prev]

# Save data
save(bp_cases_nonr, file = "bp_cases_prs_iid_icd10_filtered_all_nonrel.RData")
save(icd_cases_nonr_prev, file = "icd10_iid_cases_cbp_prev_filtered_all_nonrel.RData")
save(prev_icd_cases_f, file = "prev_of_icd10_cases_all_nonrel_filtered.RData")
save(opcs_cases_nonr_prev, file = "opcs_iid_cases_cbp_prev_filtered_all_nonrel.RData")
save(prev_opcs_cases_f, file = "prev_opcs_cases_all_nonrel_filtered.RData")
