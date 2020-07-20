#!/bin/bash
# Aim of this script is to unify original pan traits from SA Replication study using GWAS-MAP

for pain in 'Hip' 'Back' 'Neck' 'Knee' 'Head' 'Stom'

do

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/mv_gwas/data/chronic_replication/SA_10102018/scaled_filtered/MV_${pain}_Repl.SA_gwas.BGEN.stats.txt \
  --mapping-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/SA/00_original_traits/sa_repl_mapping.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/SA/00_original_traits/descriptors/chronic_${pain}_sa_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/SA/00_original_traits/unification_results/${pain}/ \
  --qc-report \
  --output-file=${pain}_output

done

