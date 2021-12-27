#!/bin/bash

path='/mnt/polyomica/projects/bp-sh/data/08_bp-sh_ma_func_an/smr-heidi/'

cat gtex_v7.yaml > ${path}/gtex_run.yaml
cat start_info.yaml >> ${path}/gtex_run.yaml

#conda activate smr-heidi-prod  # in /home/epakhomov/soft/anaconda3
conda activate /home/common/projects/glycomics/SMR-HEIDI/EUR_MA_20210826/yakov_env

PLINK_EXECUTABLE=plink /home/common/projects/glycomics/SMR-HEIDI/smr-heidi/main.py --concurrency 3 \
--create-ld-data /home/common/DataStorage/1kg_db_phase3_plink/ ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_nodup \
--skip-existing-ld-data \
-P 1 \
-r2 0.9 \
--no-png \
${path}/snps_to_run.txt ${path}/gtex_run.yaml ${path}/ld_matrices_mga ${path}/gtex_v7

conda deactivate
