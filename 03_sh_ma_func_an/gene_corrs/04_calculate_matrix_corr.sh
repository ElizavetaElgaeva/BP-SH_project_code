#to run this script it is necessary to set a filename with gwas_ids (should be listed in column) as the command line argument

array=($(grep -Pe '\d+' $1))
for gwas_id1 in ${!array[@]}
	do
		for (( gwas_id2=gwas_id1+1; gwas_id2<${#array[@]}; gwas_id2++ ))
			do
				run_ldscore \
				--rg \
				--gwas-id-1 ${array[gwas_id1]} \
				--gwas-id-2 ${array[gwas_id2]}
			done
	done 
 
