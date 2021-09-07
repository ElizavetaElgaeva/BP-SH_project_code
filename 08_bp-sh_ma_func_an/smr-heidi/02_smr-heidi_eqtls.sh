# Aim of this script is to run SMR-HEIDI analysis for BP-SH and number of eQTL traits listed in a file

while IFS= read -r line; do

	        echo "Text read from file: $line"
		run_smrheidi \
			--gwas-id-1 2269987 \
			--gwas-id-2 $line \
			--snp-list /mnt/polyomica/projects/bp-sh/data/08_bp-sh_ma_func_an/smr-heidi/snp_list.csv \
			--version 1

done < /mnt/polyomica/projects/bp-sh/data/08_bp-sh_ma_func_an/smr-heidi/locus_chr10_sh_overlap.txt

