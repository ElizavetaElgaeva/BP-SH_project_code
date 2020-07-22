#!/bin/bash
# Aim of this script is to unify Shared Heredity trait from SA Replication study using GWAS-MAP

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/bp-sh/data/01_sh/sh_replication/SA/00_raw_data/SH_SA_repl_GWAS.txt \
  --mapping-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/01_sh/SH_mapping.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/SA/01_sh/SH_SA_repl_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/SA/01_sh/unification_results/ \
  --qc-report \
  --output-file=SH_sa_repl_output

