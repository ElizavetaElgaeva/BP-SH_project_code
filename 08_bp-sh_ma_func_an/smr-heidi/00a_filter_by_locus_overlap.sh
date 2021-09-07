# Aim of this script is to filter list of eQTL traits for SMR-HEIDI analysis
# by keeeping only those which have overlapping locus with replicated locus in BP-SH

touch $1

while IFS= read -r line; do

	echo "Text read from file: $line"
	clickhouse-client -h 172.25.8.65 --query="select gwas_id from geliphe.snp where chr==10 and bp>73798873-250000 and bp<73798873+250000 and gwas_id==$line;" >> $1

done < /mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi/eqtls_ids.txt


