# Aim of this script is to recode PRS as bnary traits

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

load("./bp_prs_iid_icd10_filtered.RData")

# Recode BP PRS: 1 if the value exceeds 0.75% quantile, otherwise 0
bp_f$prs_bp_bin <- as.integer(as.logical(bp_f$"prs_bp" > quantile(bp_f$"prs_bp")[4])) 

table(bp_f$"prs_bp" > quantile(bp_f$"prs_bp")[4])
table(bp_f$prs_bp_bin)

# Recode SH PRS: 1 if the value exceeds 0.75% quantile, otherwise 0
bp_f$prs_sh_bin <- as.integer(as.logical(bp_f$"prs_sh" > quantile(bp_f$"prs_sh")[4]))

table(bp_f$"prs_sh" > quantile(bp_f$"prs_sh")[4])
table(bp_f$prs_sh_bin)

# Recode BP-SH PRS: 1 if the value exceeds 0.75% quantile, otherwise 0
bp_f$prs_bp_sh_bin <- as.integer(as.logical(bp_f$"prs_bp_sh" > quantile(bp_f$"prs_bp_sh")[4]))

table(bp_f$"prs_bp_sh" > quantile(bp_f$"prs_bp_sh")[4])
table(bp_f$prs_bp_sh_bin)

bp_f_bin <- bp_f
rm(bp_f)
reafme_bp_f_bin <- "bp_prs_iid_icd10_filtered.RData table with BP/SH/BP-SH PRS recoded as binary traits based on being into the 4th quantile"
save(bp_f_bin, reafme_bp_f_bin, file = "./bp_prs_bin_iid_filtered.RData")

