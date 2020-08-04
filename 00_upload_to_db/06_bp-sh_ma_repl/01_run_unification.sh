#!/bin/bash
# Aim of this script is to unify Replication meta-analysis of Back pain without Shared Heredity trait using GWAS-MAP

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/bp-sh/data/05_bp-sh_ma/02_repl_ma/02_ma_METAL/BP-SH_ma_repl_metal.txt \
  --mapping-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/05_bp-sh_ma_eur/BP-SH_ma_mapping.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/06_bp-sh_ma_repl/BP-SH_ma_repl_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/06_bp-sh_ma_repl/unification_results/ \
  --qc-report \
  --output-file=BP-SH_ma_repl_output

