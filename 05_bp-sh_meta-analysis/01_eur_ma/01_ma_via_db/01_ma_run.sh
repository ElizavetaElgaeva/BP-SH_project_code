# Aim of this script is to run EUR MA of Back pain without Shared heredity using GWAS-MAP

out="/mnt/polyomica/projects/bp-sh/data/05_bp-sh_ma/01_eur_ma/01_ma_via_db/"

run_meta_analysis --gwas-ids 41,42 --genomic-control ON --output-folder "$out" --output-filename BP-SH_EUR_MA.csv


