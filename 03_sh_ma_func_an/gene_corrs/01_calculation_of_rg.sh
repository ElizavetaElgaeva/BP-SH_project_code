# This script is for calculation of genetic correlations between EUR_MA/EUR_MA-SH and all of the non eQTL GWASes with sample size>10000 (h2>0.001 and Z-score(h2)>3)

traits_of_interest='2269987'
gwas_id=$(cat ../../../../data/03_sh_ma_func_an/gene_corrs/20200911_list_of_gwas_ids_for_rg.txt)
for trait in $traits_of_interest
do
	echo Calculations started for trait id $trait
	for gid in $gwas_id
	do
		#echo $trait $gid
		run_ldscore \
		--rg \
		--gwas-id-1 $trait \
		--gwas-id-2 $gid
	done
done
 
