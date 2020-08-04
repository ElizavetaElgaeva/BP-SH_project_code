#!/bin/bash
# Aim of this script is to unify European meta-analysis of Back pain without Shared Heredity trait using GWAS-MAP

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/bp-sh/data/05_bp-sh_ma/01_eur_ma/02_ma_METAL/BP-SH_ma_eur_metal.txt \
  --mapping-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/05_bp-sh_ma_eur/BP-SH_ma_mapping.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/05_bp-sh_ma_eur/BP-SH_ma_eur_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/05_bp-sh_ma_eur/unification_results/ \
  --qc-report \
  --output-file=BP-SH_ma_eur_output

