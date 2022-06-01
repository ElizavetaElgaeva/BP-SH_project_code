# Aim of this script is to filter OPCS data

library(data.table)
library(readxl)
library(xlsx)

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

# Recalculate phenotypes using combined OPCS codes 
opcs_full <- read.xlsx("OPCS_coding240.xlsx", sheetIndex = 1) # read full OPCS code table, definition and description of each code 

level_2 <- sapply(opcs_full$coding, function(x) substr(x, 1, 3)) # extract first 3 symbols from coicd_full

level_2 <- unique(level_2)
length(level_2) # 1524

codes_to_combine <- sapply(level_2, function(x) grep(colnames(opcs), pattern = x)) # select which OPCS codes from UKBB data could be combined

grouped_pheno <- lapply(codes_to_combine, function(x) opcs_f[, x]) # list of grouped phenotypes by OPCS codes of the 2nd level
grouped_pheno <- lapply(grouped_pheno, function(x) as.data.frame(x)) 

pulled_pheno <- lapply(grouped_pheno, function(x) rowSums(x)) # obtain the sum of all columns with codes within one level
pulled_pheno <- lapply(pulled_pheno, function(x) replace(x, !is.na(x) & x != 0, 1)) # recode the sum in binary way

tmp <- as.data.frame(do.call(cbind, pulled_pheno)) # obtain a data frame from list of sum vectors
opcs_f <- opcs_f[, c(1, 8003:8023)] # extract part of the table containing no OPCS codes
opcs_f_l2 <- cbind(opcs_f, tmp)

# Save data
readme_opcs_f_l2 <- "OPCS table (level 2 codes) with only IIDs presented in bp_prs_iid_icd10_filtered.RData"
save(opcs_f_l2, readme_opcs_f_l2, file = "./opcs_level_2_iid_cbp_filtered.RData")


