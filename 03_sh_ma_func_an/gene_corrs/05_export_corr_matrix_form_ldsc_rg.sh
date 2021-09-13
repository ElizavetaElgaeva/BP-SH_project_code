# Export correlation matrix from ldsc_rg table\
# Before start this script it is necessary to calculate rg for the input
# file containing GWAS ids.
# To start this scirpt it is necessary to set path to the above 
# mentioned file as the first argument of command line
# the second argument should contain path to output file
#query=''
array=($(grep -Pe '\d+' $1))

#Query to extract headers without any data
test_psql_connect -c "COPY (SELECT * FROM gwas.ldsc_rg WHERE false) TO STDOUT With CSV DELIMITER ',' HEADER" > $2
for gwas_id1 in ${!array[@]}
	do
		for (( gwas_id2=gwas_id1+1; gwas_id2<${#array[@]}; gwas_id2++ ))
			do
				prod_psql_connect -c "COPY (SELECT * FROM gwas.ldsc_rg WHERE gwas_id_1=${array[gwas_id1]} AND gwas_id_2=${array[gwas_id2]}) TO STDOUT With CSV DELIMITER ','" >> $2
			done
	done 
