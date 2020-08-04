#!/bin/bash
# Aim of this script is to unify European meta-analysis of GIP1
# for four musculosceletal traits from our previous study (https://doi.org/10.1038/s42003-020-1051-9) using GWAS-MAP

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/GPC1/gpc1_Pval \
  --mapping-path=/mnt/polyomica/projects/mv_gwas/MA/GPC_chronic_MA_eur/upload_to_db/mapping_gpc_ma_chron_pain_eur.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/02_gip1_eur_ma_previous_study/gpc1_chron_ma_eur_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/02_gip1_eur_ma_previous_study/unification_results/ \
  --qc-report \
  --output-file=gip1_ma_eur_output

