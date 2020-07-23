# Aim of this script is to run EUR MA of Shared heredity using GWAS-MAP

out="/mnt/polyomica/projects/bp-sh/data/02_sh_ma/01_eur_ma/01_ma_via_db/"

run_meta_analysis --gwas-ids 19,26 --genomic-control ON --output-folder "$out" --output-filename SH_EUR_MA.csv


