#!/bin/bash
# Aim of this script is to run COJO on SH MA Eur data with p = 5e-06 threshold

path='/storage/projects/PainOmics/bp-sh/data/cojo/01_sh/'

for i in {1..22}
do
   /storage/projects/PainOmics/BP_crude_GWAS/gcta_1.90.0beta/gcta64 \
   --bfile /storage/projects/PainOmics/MV_GWAS/100k_bed_filtered/100k_chr"$i" \
   --maf 0.0002 \
   --cojo-slct \
   --cojo-p 5e-6 \
   --chr $i \
   --cojo-wind 5000 \
   --cojo-file ${path}/input/SH_ma_eur_output_done.for_cojo \
   --out ${path}/output/5e-6/sh_ma_eur_chr"$i"_5e-6.cojo
done
