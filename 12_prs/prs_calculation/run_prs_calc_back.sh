#!/bin/bash
path_to_gwas="/home/anostaeva/Back_Pain_2021/gwas/Back_3000_genome.csv"
plink2="/home/anostaeva/plink2/plink2"
path_to_list="/home/anostaeva/genotype/ukb_1000145/pgen/stat"

for i in {1..22}
do
plink2 --pfile /home/anostaeva/genotype/ukb_1000145/pgen/ukb_imp_chr${i}_v3_filtered \
       --keep ${path_to_list}/individuals_filtered_0.98.txt \
       --extract ${path_to_list}/snps_filtered_0.98.txt \
       --score ${path_to_gwas} 1 3 6 header cols=scoresums \
       --out /home/anostaeva/Back_Pain_2021/UKB_all_PRS/Back/ukb_imp_chr_${i}_v3_PRS
done

python3 create_prs.py

rm /home/anostaeva/Back_Pain_2021/UKB_all_PRS/Back/*sscore*
rm /home/anostaeva/Back_Pain_2021/UKB_all_PRS/Back/*log*

