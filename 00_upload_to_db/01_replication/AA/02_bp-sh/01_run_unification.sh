#!/bin/bash
# Aim of this script is to unify Back pain without Shared Heredity trait from AA Replication study using GWAS-MAP

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/bp-sh/data/04_bp-sh/bp-sh_replication/AA/BP-SH_aa_repl.txt \
  --mapping-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/02_bp-sh/BP-SH_mapping.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/AA/02_bp-sh/BP-SH_AA_repl_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/AA/02_bp-sh/unification_results/ \
  --qc-report \
  --output-file=BP-SH_aa_repl_output

