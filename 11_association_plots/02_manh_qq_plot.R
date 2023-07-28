# Aim of this script is to build Manhattan and QQ plots
# for SGIT and CBP UGIT (discovery sample)

library(data.table)
library(qqman)
library(qqman)

setwd("/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/")

path_to_fig <- "/mnt/polyomica/projects/bp-sh/data/11_association_plots/"

file_names <- c("./01_sh/unification_results/SH_disc_output_done.csv", "./02_bp-sh/unification_results/BP-SH_disc_output_done.csv")

intercept <- c(1.03758691779465, 0.987876106164573) # LD Score intercept for SGIT and UGIT Discovery calculated using --h2 option through the data base

for(i in 1:lenght(file_names)){

tmp <- fread(file_names[i])

gwas_name <- tail(strsplit(file_names[i], split = "/")[[1]], n = 1)
gwas_name <- gsub("_output_done.csv", "", gwas_name)
short_name <- gsub("_disc", "", gwas_name)

# Estimate lambda GC
z2 <- (tmp$z)^2
lambda <- median(z2)/qchisq(0.5,1)
lambda <- round(lambda, 4)

# Build a Q-Q plot
tiff(paste0(path_to_fig, short_name, "/qq_plot_", gwas_name, ".tiff"))
qq(tmp$p, main = paste0("lambda = ", lambda))
dev.off()

# Correction for genomic control
tmp$p_gc <- pchisq(((tmp$beta / tmp$se)^2 / intercept[i]), df = 1, lower.tail = F)

# Filter by MAF
tmp$maf <- pmin(tmp$eaf, 1 - tmp$eaf)
tmp <- tmp[tmp$maf >= 2e-04, ]

# Creae a Manhattan plot
tiff(paste0(path_to_fig, short_name, "/manh_plot_", gwas_name,".tiff"), height = 1080, width = 1920, res = 150, pointsize = 8)
manhattan(tmp, chr="chr", bp="bp", snp="rs_id", p="p_gc", genomewideline = -log10(5e-08/6), suggestiveline = F, ylim = c(0, 15), col = c("red2", "steelblue3", "green3", "mediumorchid3", "darkorange1", "yellow1", "chocolate", "hotpink1", "lavenderblush3"))
dev.off()

}
