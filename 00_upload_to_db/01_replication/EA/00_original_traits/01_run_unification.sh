#!/bin/bash
# Aim of this script is to unify original pan traits from EA Replication study using GWAS-MAP

for pain in 'Hip' 'Back' 'Neck' 'Knee' 'Head' 'Stom'

do

run_uni_qc_rep \
  --gwas-path=/mnt/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/scaled_filterted/MV_${pain}_Repl.EA_gwas.BGEN.stats.txt \
  --mapping-path=/mnt/polyomica/projects/mv_gwas/data/chronic_replication/EA_09102018/upload_to_db/mapping_chronic_pain_replication.json \
  --descriptors-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/EA/00_original_traits/descriptors/chronic_${pain}_eur_descriptor.json \
  --output-dir=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/EA/00_original_traits/unification_results/${pain}/ \
  --qc-report \
  --output-file=${pain}_output

done

