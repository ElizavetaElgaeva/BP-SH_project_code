# Aim of this script is to run SMR-HEIDI analysis for SH and number of comlex traits listed in a file

for num in 1 3 4 7 9 13 17; do

while IFS= read -r line; do

	        echo "Text read from file: $line"
		run_smrheidi \
			--gwas-id-1 2269986 \
			--gwas-id-2 $line \
			--snp-list /mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi/chr${num}.csv \
			--version 1

done < /mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi/locus_chr${num}_sh_overlap.txt

done

