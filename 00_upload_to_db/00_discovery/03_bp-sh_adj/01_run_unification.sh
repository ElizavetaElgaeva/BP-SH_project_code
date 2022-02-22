#!/bin/bash
# Aim of this script is to unify adjusted CBP UGIT from Discovery study using GWAS-MAP

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/bp-sh/data/13_bp-sh_adj/bp-sh_adj_discovery/BP-SH_adj_discovery.txt \
  --mapping-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/03_bp-sh_adj/BP-SH_adj_mapping.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/03_bp-sh_adj/BP-SH_adj_discovery_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/03_bp-sh_adj/unification_results/ \
  --qc-report \
  --output-file=BP-SH_adj_disc_output

