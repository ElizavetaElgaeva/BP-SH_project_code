# Aim of this script is to standardize PRS values

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS")

# Load PRS data (test sample, non-relatieves only)
load("bp_prs_iid_icd10_filtered_test_nonrelatives.RData")
ls()

colnames(bp_f_test_nonr)
bp_f_test_nonr$"prs_bp" <- bp_f_test_nonr$"prs_bp"/sqrt(var(bp_f_test_nonr$"prs_bp"))
bp_f_test_nonr$"prs_sh" <- bp_f_test_nonr$"prs_sh"/sqrt(var(bp_f_test_nonr$"prs_sh"))
bp_f_test_nonr$"prs_bp_sh" <- bp_f_test_nonr$"prs_bp_sh"/sqrt(var(bp_f_test_nonr$"prs_bp_sh"))

bp_f_test_nonr_var1 <- bp_f_test_nonr

save(bp_f_test_nonr_var1, file = "bp_prs_iid_icd10_filtered_test_nonrelatives_prs_var1.RData")
