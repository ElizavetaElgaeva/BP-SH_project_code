#!/bin/bash
# Aim of this script is to unify original pain traits from Discovery study using GWAS-MAP

for pain in 'Hip' 'Back' 'Neck' 'Knee' 'Head' 'Stom'

do

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/scaled_filtered/MV_${pain}_Disc_gwas.BGEN.stats.txt \
  --mapping-path=/mnt/polyomica/projects/mv_gwas/data/chronic_discovery/upload_to_db/mapping_chronic_pain_discovery.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/00_original_traits/descriptors/chronic_${pain}_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/00_discovery/00_original_traits/unification_results/${pain}/ \
  --qc-report \
  --output-file=${pain}_output

done

