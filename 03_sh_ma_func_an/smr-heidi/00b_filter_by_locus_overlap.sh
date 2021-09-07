# Aim of this script is to filter list of eQTL traits for SMR-HEIDI analysis
# by keeeping only those which have overlapping locus with replicated locus in SH

touch $3

while IFS= read -r line; do

	echo "Text read from file: $line"
	clickhouse-client -h 172.25.8.65 --query="select gwas_id from geliphe.snp where chr==$1 and bp>$2-250000 and bp<$2+250000 and gwas_id==$line;" >> $3

done < /mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi/eqtls_ids.txt


