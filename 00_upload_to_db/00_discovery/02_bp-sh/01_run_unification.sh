#!/bin/bash
# Aim of this script is to unify Back pain without Shared Heredity trait from Discovery study using GWAS-MAP

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/bp-sh/data/04_bp-sh/bp-sh_discovery/BP-SH_discovery.txt \
  --mapping-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/02_bp-sh/BP-SH_mapping.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/02_bp-sh/BP-SH_discovery_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/02_bp-sh/unification_results/ \
  --qc-report \
  --output-file=BP-SH_disc_output

