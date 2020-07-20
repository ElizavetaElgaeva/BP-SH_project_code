#!/bin/bash
# Aim of this script is to upload original pain traits from AA Replication study to GWAS-MAP

for pain in 'Hip' 'Back' 'Neck' 'Knee' 'Head' 'Stom'

do

run_upload --gwas-path=/mnt/polyomica/projects/bp-sh/data/00_upload_to_db/01_replication/AA/00_original_traits/unification_results/${pain}/${pain}_output_done.csv

done

