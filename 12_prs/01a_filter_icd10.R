# Aim of this script is to filter ICD10 data

library(data.table)

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/")

# Read table with ukbb data ICD10 codes 
icd <- fread("./ICD10_main_secondary_merged.txt", data.table = F)

# Read clean table for BP
bp <- fread("./cbp_pheno_prs_merged.txt", data.table = F)

# Filter ICD10 data to keep only those IIDs, that are in clean BP table
nrow(icd)
#[1] 487322
nrow(bp)
#[1] 439829
table(bp$IID %in% icd$IID)
# FALSE   TRUE
#    67 439762
table(icd$IID %in% bp$IID)
# FALSE   TRUE
# 47560 439762

ind <- match(bp$IID, icd$IID)
ind <- na.omit(ind)
icd_f <- icd[ind, ]

table(icd_f$IID %in% bp$IID)

ind2 <- match(icd_f$IID, bp$IID)
bp_f <- bp[ind2, ] # filter CBP, PRS table as well
table(icd_f$IID == bp_f$IID)

# Check if there non-binary traits
# columns 1 and the last 21 are not ICD10 codes
bin <- lapply(icd_f, function(x) { all(x %in% 0:1) })
table(unlist(bin)[-c(1, 11471:11491)] == TRUE)
# TRUE
#11469
# no non-binary traits

# Check for NAs
na <- lapply(icd_f, function(x) all(!is.na(x)))
table(unlist(na) == TRUE)
#FALSE  TRUE
#    7 11484
which(unlist(na) == FALSE)
#              Headache             Facialpain     Neckorshoulderpain
#                 11472                  11473                  11474
#              Backpain Stomachorabdominalpain                Hippain
#	                       11475                  11476                  11477
#	                    Kneepain
#			                     11478

dim(na.omit(icd_f))
#[1] 437604  11491

# Check if empty
empty <- lapply(icd_f, function(x) all(x != ""))
table(unlist(empty) == TRUE)
# TRUE
#11484
# no emty fields

# Filter by prevalence
prev <- lapply(icd_f[, -c(1, 11471:11491)], function(x) sum(x))
prev <- unlist(prev)
n <- nrow(icd_f)
prev <- lapply(prev, function(x) x/n*100)
ind_prev <- which(prev > 0.5 & prev < 99.5)
icd_f_prev <- icd_f[, c(1, ind_prev + 1)]

# Save data
icd_f <- icd_f[, -c(11471:11491)]
readme_icd_f <- "ICD10 table with only IIDs presented in cbp_pheno_prs_merged.txt"
save(icd_f, readme_icd_f, file = "./icd10_iid_cbp_filtered.RData")

readme_bp_f <- "CBP PRS table with only IIDs presented in icd10_iid_cbp_filtered.RData"
save(bp_f, readme_bp_f, file = "./bp_prs_iid_icd10_filtered.RData")

readme_icd_f_prev <- "ICD10 table with only IIDs presented in cbp_pheno_prs_merged.txt and ICD10 codes for binary traits with prevalence > 0.5% and < 99.5%"
save(icd_f_prev, readme_icd_f_prev, file = "./icd10_iid_cbp_prev_filtered.RData")


