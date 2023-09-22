#!/bin/bash
path_to_gwas="/home/anostaeva/Back_Pain_project/Back_Pain_2022/gwas"
path_to_out="/home/anostaeva/Back_Pain_project/Back_Pain_2022/UKB_all_PRS/bp_sh"
plink2="/home/anostaeva/plink2/plink2"
path_to_list="/home/anostaeva/genotype/ukb_1000145/stat"
path_to_ukb="/home/anostaeva/genotype/ukb_1000145/pgen"
path_to_code="/home/anostaeva/Back_Pain_project/Back_Pain_2022/UKB_all_PRS/bp_sh"

cd /home/anostaeva

for id in {1..6}; do
for i in {1..22}; do

plink2 --pfile ${path_to_ukb}/ukb_imp_chr${i}_v3_filtered \
       --keep ${path_to_list}/individuals_filtered_0.98.txt \
       --extract ${path_to_list}/snps_filtered_0.98.txt \
       --score ${path_to_gwas}/BP_SH_model_${id}_genome.csv 1 3 6 header cols=scoresums \
       --out ${path_to_out}/ukb_imp_chr_${i}_v3_PRS
done

singularity exec /home/anostaeva/py.sif python3 ${path_to_code}/create_prs.py ${path_to_out} ukb_bp_sh_${id}

rm ${path_to_out}/*sscore*
rm ${path_to_out}/*log*

done
