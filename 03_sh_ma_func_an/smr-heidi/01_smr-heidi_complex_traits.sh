# Aim of this script is to run SMR-HEIDI analysis for SH and number of comlex traits listed in a file

while IFS= read -r line; do

	echo "Text read from file: $line"
	run_smrheidi \
		--gwas-id-1 2269986 \
		--gwas-id-2 $line \
		--snp-list /mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi/snp_list.csv \
		--version 1

done < /mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi/complex_traits_ids.txt


