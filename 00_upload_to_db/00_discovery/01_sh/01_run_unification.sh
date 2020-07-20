#!/bin/bash
# Aim of this script is to unify Shared Heredity trait from Discovery study using GWAS-MAP

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/bp-sh/data/01_sh/sh_discovery/00_raw_data/SH_disc_GWAS.txt \
  --mapping-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/01_sh/SH_mapping.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/01_sh/SH_discovery_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/01_sh/unification_results/ \
  --qc-report \
  --output-file=SH_disc_output

