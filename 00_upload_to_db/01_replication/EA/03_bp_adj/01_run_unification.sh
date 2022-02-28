#!/bin/bash
# Aim of this script is to unify adjusted CBP from EA Replication study using GWAS-MAP

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/bp-sh/data/13_bp_adj/bp_adj_replication/EA/BP_adj_ea_repl.txt \
  --mapping-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/02_bp-sh/BP-SH_mapping.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/EA/03_bp_adj/BP_adj_EA_repl_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/EA/03_bp_adj/unification_results/ \
  --qc-report \
  --output-file=BP_adj_ea_repl_output

