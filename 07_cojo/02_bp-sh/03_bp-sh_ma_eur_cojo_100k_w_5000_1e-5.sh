#!/bin/bash
# Aim of this script is to run COJO on BP-SH MA Eur data with p = 1e-05 threshold

path='/storage/projects/PainOmics/bp-sh/data/cojo/02_bp-sh/'

for i in {1..22}
do
   /storage/projects/PainOmics/BP_crude_GWAS/gcta_1.90.0beta/gcta64 \
   --bfile /storage/projects/PainOmics/MV_GWAS/100k_bed_filtered/100k_chr"$i" \
   --maf 0.0002 \
   --cojo-slct \
   --cojo-p 1e-5 \
   --chr $i \
   --cojo-wind 5000 \
   --cojo-file ${path}/input/BP-SH_ma_eur_output_done.for_cojo \
   --out ${path}/output/1e-5/bp-sh_ma_eur_chr"$i"_1e-5.cojo
done
