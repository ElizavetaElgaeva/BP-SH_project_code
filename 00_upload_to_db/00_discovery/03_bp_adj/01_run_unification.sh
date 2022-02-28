#!/bin/bash
# Aim of this script is to unify adjusted CBP from Discovery study using GWAS-MAP

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/bp-sh/data/13_bp_adj/bp_adj_discovery/BP_adj_discovery.txt \
  --mapping-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/03_bp_adj/BP_adj_mapping.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/03_bp_adj/BP_adj_discovery_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/03_bp_adj/unification_results/ \
  --qc-report \
  --output-file=BP_adj_disc_output

