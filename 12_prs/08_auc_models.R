# Aim of this script is to estimate AUC value for models using PRS for BP/SH/BP-SH
# Test sample, non-relatives only; standardized PRS

library(ModelMetrics)

setwd("/home/common/projects/pain_project/UKBB_pheno_ICD10_OPCS/")

# Load PRS data
load("bp_prs_iid_icd10_filtered_test_nonrelatives_prs_var1.RData")

m1 <- glm(back ~ prs_bp + Age + Sex, data = bp_f_test_nonr_var1, family = "binomial")
auc(m1)
# 0.5625393

m2 <- glm(back ~ prs_bp + Age + Sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, data = bp_f_test_nonr_var1, family = "binomial")
auc(m2)
# 0.5635858

m3 <- glm(back ~ prs_bp*bmi + Age + Sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, data = bp_f_test_nonr_var1, family = "binomial")
auc(m3)
# 0.5900467

m4 <- glm(back ~ prs_bp*prs_sh*prs_bp_sh*bmi + Age + Sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, data = bp_f_test_nonr_var1, family = "binomial")
auc(m4)
# 0.5947017

m5 <- glm(back ~ prs_sh + Age + Sex, data = bp_f_test_nonr_var1, family = "binomial")
auc(m5)
# 0.5696085

m6 <- glm(back ~ prs_sh + Age + Sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, data = bp_f_test_nonr_var1, family = "binomial")
auc(m6)
# 0.5702607

m7 <- glm(back ~ prs_sh*bmi + Age + Sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, data = bp_f_test_nonr_var1, family = "binomial")
auc(m7)
# 0.5929865

m8 <- glm(back ~ prs_bp_sh + Age + Sex, data = bp_f_test_nonr_var1, family = "binomial")
auc(m8)
# 0.518918

m9 <- glm(back ~ prs_bp_sh + Age + Sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, data = bp_f_test_nonr_var1, family = "binomial")
auc(m9)
# 0.5223091

m10 <- glm(back ~ prs_bp_sh*bmi + Age + Sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, data = bp_f_test_nonr_var1, family = "binomial")
auc(m10)
# 0.5709563




