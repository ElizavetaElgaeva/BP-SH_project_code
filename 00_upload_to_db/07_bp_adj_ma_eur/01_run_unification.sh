#!/bin/bash
# Aim of this script is to unify European meta-analysis of Back pain adjusted for all other original pain traits using GWAS-MAP

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/bp-sh/data/14_bp_adj_ma/01_eur_ma/01_ma_METAL/BP_adj_ma_eur_metal.txt \
  --mapping-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/05_bp-sh_ma_eur/BP-SH_ma_mapping.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/07_bp_adj_ma_eur/BP_adj_ma_eur_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/07_bp_adj_ma_eur/unification_results/ \
  --qc-report \
  --output-file=BP_adj_ma_eur_output

