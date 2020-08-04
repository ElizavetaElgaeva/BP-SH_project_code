#!/bin/bash
# Aim of this script is to unify European meta-analysis of Shared Heredity trait using GWAS-MAP

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/bp-sh/data/02_sh_ma/01_eur_ma/02_ma_METAL/SH_ma_eur_metal.txt \
  --mapping-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/03_sh_ma_eur/SH_ma_mapping.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/03_sh_ma_eur/SH_ma_eur_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/03_sh_ma_eur/unification_results/ \
  --qc-report \
  --output-file=SH_ma_eur_output

