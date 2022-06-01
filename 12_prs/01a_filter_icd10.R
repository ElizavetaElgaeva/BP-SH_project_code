# Aim of this script is to filter ICD10 data

library(data.table)
library(readxl)
library(xlsx)

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

# Recalculate phenotypes using combined ICD10 codes 
icd_full <- read.xlsx("ICD10_coding19.xlsx", sheetIndex = 1) # read full ICD10 code table, definition and description of each code 

level_2 <- sapply(icd_full$coding, function(x) substr(x, 1, 3)) # extract first 3 symbols from code

level_2 <- unique(level_2)
length(level_2) # 2062

codes_to_combine <- sapply(level_2, function(x) grep(colnames(icd), pattern = x)) # select which ICD10 codes from UKBB data could be combined

grouped_pheno <- lapply(codes_to_combine, function(x) icd_f[, x]) # list of grouped phenotypes by ICD10 codes of the 2nd level
grouped_pheno <- lapply(grouped_pheno, function(x) as.data.frame(x)) 

pulled_pheno <- lapply(grouped_pheno, function(x) rowSums(x)) # obtain the sum of all columns with codes within one level
pulled_pheno <- lapply(pulled_pheno, function(x) replace(x, !is.na(x) & x != 0, 1)) # recode the sum in binary way

tmp <- as.data.frame(do.call(cbind, pulled_pheno)) # obtain a data frame from list of sum vectors
icd_f <- icd_f[, c(1, 11471:11491)]
icd_f_l2 <- cbind(icd_f, tmp)

# Save data
readme_icd_f_l2 <- "UKBB ICD10 table with codes combined to level 2. Not filtered by prevalence. Full sample size, N=439K. Only IIDs presented in cbp_pheno_prs_merged.txt"
save(icd_f_l2, readme_icd_f_l2, file = "./icd10_iid_level_2_cbp_filtered.RData")

readme_bp_f <- "CBP PRS table with only IIDs presented in icd10_iid_cbp_filtered.RData"
save(bp_f, readme_bp_f, file = "./bp_prs_iid_icd10_filtered.RData")


