# Aim of this script is to filter the results of SMR-HEIDI analysis for BP-SH

library(data.table)

setwd("/mnt/polyomica/projects/bp-sh/data/08_bp-sh_ma_func_an/smr-heidi")

# Read a file with preselected tissue types that are needed
tissues <- fread("/mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi/matching_tissues.txt", data.table = F)

# Read the SMR-HEIDI results
gtex <- fread("gtex_v7/results.tsv", data.table = F)
dim(gtex)
cedar <- fread("cedar/results.tsv", data.table = F)
dim(cedar)
westra <- fread("westra/results.tsv", data.table = F)

# A function to extract certain tissue types from SMR-HEIDI reults using predefined list of tissues
filter_smar_heidi_results <- function(tissues, smr_heidi_res){
	matches <- lapply(tissues, function(x) grep(x, smr_heidi_res$gwas_2))
	matches <- unlist(matches)
	filterted <- smr_heidi_res[matches, ]
	return(filterted)
}

# Filter GTEX data
gtex_f <- filter_smar_heidi_results(tissues = tissues$tissue_name[tissues$Dataset == "GTEX_v7" & !is.na(tissues$tissue_name)], smr_heidi_res = gtex) 
gtex_f$Source <- "GTEx_v7"

# Filter CEDAR data
cedar_f <- filter_smar_heidi_results(tissues = tissues$tissue_name[tissues$Dataset == "cedar" & !is.na(tissues$tissue_name)], smr_heidi_res = cedar)
cedar_f$Source <- "CEDAR"

# Westra should not be filtered
westra$Source <- "Westra"

# Merge and save
# Check columns
length(colnames(gtex_f)) # 53
length(colnames(cedar_f)) # 52
length(colnames(westra)) # 53
which(!(colnames(gtex_f) %in% colnames(cedar_f))) # 46
colnames(gtex_f)[46] # maf_G2
table(colnames(gtex_f)[-46] == colnames(cedar_f)) # all true
table(colnames(gtex_f) == colnames(westra)) # all true
# can ignor maf_G2 column
merged_full <- rbind(gtex_f[, -c(46)], cedar_f)
merged_full <- rbind(merged_full, westra[, -c(46)])
fwrite(merged_full, file = "bp-sh_results_full_preselected_tissues.txt", dec = ".", sep = "\t")

# Let's filter the columns
merged_filtered <- merged_full[, c("locus", "rsid_REF", "locus_top_ld_r", "gwas_1", "gwas_2", "Source", "gene_name_G2", "beta_smr", "sigma_smr", "p_smr", "p_heidi", "n_heidi", "p_G2", "probe_G2")]
fwrite(merged_filtered, file = "bp-sh_results_filtered_columns.txt", dec = ".", sep = "\t")

# There is no replicated SNP
