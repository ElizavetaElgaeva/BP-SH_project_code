#!/bin/bash

path_to_phenotype="/home/anostaeva/Back_Pain_project/Back_Pain_2022/phenotype/phenotypic_data.txt"
path_to_prs="/home/anostaeva/Back_Pain_project/Back_Pain_2022/UKB_all_PRS"

for fullfile in ${path_to_prs}/back/ukb_back_2.csv; do

trait="back"
echo ${trait}
echo ${fullfile}
python3 accuracy_test.py ${fullfile} ${path_to_phenotype} ${trait}

done

for fullfile in ${path_to_prs}/bp_sh/ukb_bp_sh_1.csv; do

trait="bp_sh"
echo ${trait}
echo ${fullfile}
python3 accuracy_test.py ${fullfile} ${path_to_phenotype} ${trait}

done

for fullfile in ${path_to_prs}/bp_adj/ukb_bp_adj_5.csv; do

trait="bp_adj"
echo ${trait}
echo ${fullfile}
python3 accuracy_test.py ${fullfile} ${path_to_phenotype} ${trait}

done

for fullfile in ${path_to_prs}/sh/ukb_sh_2.csv; do

trait="sh"
echo ${trait}
echo ${fullfile}
python3 accuracy_test.py ${fullfile} ${path_to_phenotype} ${trait}

done
