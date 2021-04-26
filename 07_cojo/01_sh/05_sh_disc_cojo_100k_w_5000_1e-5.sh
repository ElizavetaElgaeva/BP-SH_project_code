#!/bin/bash
# Aim of this script is to run COJO on SH Discovery data with p = 1e-05 threshold

path='/storage/projects/PainOmics/bp-sh/data/cojo/01_sh/'

for i in {1..22}
do
   /storage/projects/PainOmics/BP_crude_GWAS/gcta_1.90.0beta/gcta64 \
   --bfile /storage/projects/PainOmics/MV_GWAS/100k_bed_filtered/100k_chr"$i" \
   --maf 0.0002 \
   --cojo-slct \
   --cojo-p 1e-5 \
   --chr $i \
   --cojo-wind 5000 \
   --cojo-file ${path}/input/SH_disc_output_done.for_cojo \
   --out ${path}/output/discovery/1e-5/sh_disc_chr"$i"_1e-5.cojo
done
