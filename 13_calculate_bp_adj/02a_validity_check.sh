# Aim of this script is to estimate heritability of CBP
# adjusted for all other original pain traits
# using LD Score implementd in GWAS-Map (from test container). EA Replication study

# Estimate heritability for adjusted BP
run_ldscore --h2 --gwas-id=54
# 0.0126 ## se 0.0028 

# Estimate genetic correlation for SH (disc) and adjusted BP
run_ldscore --rg --gwas-id=19,54
## 0.2047 (se 0.06, p 0.0007)

# Estimate genetic correlation for SH (ea repl) and adjusted BP
run_ldscore --rg --gwas-id=26,54
## 0.1176 (se 0.073, p 0.106)

# Estimate genetic correlation for BP (disc) and adjusted BP
run_ldscore --rg --gwas-id=14,54
## 0.3838 (se 0.081, p 2.3233e-06)

# Estimate genetic correlation for BP (ea repl) and adjusted BP
run_ldscore --rg --gwas-id=21,54
## 0.4404 (se 0.066, p 2.9465e-11)

# Estimate genetic correlation for BP-SH (disc) and adjusted BP
run_ldscore --rg --gwas-id=41,54
## 0.4656 (se 0.13, p 0.0003)

# Estimate genetic correlation for BP-SH (ea repl) and adjusted BP
run_ldscore --rg --gwas-id=42,54
## 0.8557 (se 0.043, p 1.6090e-89)

# Estimate genetic correlation for Hip pain (disc) and adjusted BP
run_ldscore --rg --gwas-id=13,54
## 0.2833 (se 0.098, p 0.0037)

# Estimate genetic correlation for Hip pain (ea repl) and adjusted BP
run_ldscore --rg --gwas-id=20,54
## 0.0838 (se 0.123, p 0.4947)

# Estimate genetic correlation for Neck pain (disc) and adjusted BP
run_ldscore --rg --gwas-id=15,54
## 0.0244 (se 0.087, p 0.7795)

# Estimate genetic correlation for Neck pain (ea repl) and adjusted BP
run_ldscore --rg --gwas-id=22,54
## -0.2316 (se 0.086, p 0.0073)

# Estimate genetic correlation for Knee pain (disc) and adjusted BP
run_ldscore --rg --gwas-id=16,54
## 0.1772 (se 0.078, p 0.0233)

# Estimate genetic correlation for Knee pain (ea repl) and adjusted BP
run_ldscore --rg --gwas-id=23,54
## 0.1638 (se 0.1, p 0.1027)

# Estimate genetic correlation for Head pain (disc) and adjusted BP
run_ldscore --rg --gwas-id=17,54
## -0.1088 (se 0.082, p 0.1836)

# Estimate genetic correlation for Head pain (ea repl) and adjusted BP
run_ldscore --rg --gwas-id=24,54
## 0.0117 (se 0.089, p 0.8961)

# Estimate genetic correlation for Stomach pain (disc) and adjusted BP
run_ldscore --rg --gwas-id=18,54
## 0.1383 (se 0.114, p 0.2259)

# Estimate genetic correlation for Stomach pain (ea repl) and adjusted BP
run_ldscore --rg --gwas-id=25,54
## 0.0133 (se 0.131, p 0.9192)

# Estimate genetic correlation for adjusted BP (ea repl) and adjusted BP (disc)
run_ldscore --rg --gwas-id=53,54
## 0.7591 (se 0.188, p 5.5165e-05)


# compare with the results of 01a script
