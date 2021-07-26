# Aim of this script is to reformat the unified GWAS data
# on SH and BP-SH (Discovery) to make them suitable for locuszoom

library(data.table)

setwd('/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/')

sh <- fread('./01_sh/unification_results/SH_disc_output_done.csv')
sh <- sh[, c('chr', 'bp', 'ra', 'ea', 'p', 'beta', 'se', 'eaf')]
colnames(sh) <- c('Chromosome', 'Position', 'Ref. allele', 'Alt. allele', 'p-value', 'Effect size', 'Standard error of the effect size', 'Frequency')
sh <- sh[order(Chromosome, Position), ]
head(sh)

fwrite(sh, './01_sh/unification_results/SH_disc_output_for_locuszoom.txt', sep = '\t', dec = '.')

bpsh <- fread('./02_bp-sh/unification_results/BP-SH_disc_output_done.csv')
bpsh <- bpsh[, c('chr', 'bp', 'ra', 'ea', 'p', 'beta', 'se', 'eaf')]
colnames(bpsh) <- c('Chromosome', 'Position', 'Ref. allele', 'Alt. allele', 'p-value', 'Effect size', 'Standard error of the effect size', 'Frequency')
bpsh <- bpsh[order(Chromosome, Position), ]
head(bpsh)

fwrite(bpsh, './02_bp-sh/unification_results/BP-SH_disc_output_for_locuszoom.txt', sep = '\t', dec = '.')

