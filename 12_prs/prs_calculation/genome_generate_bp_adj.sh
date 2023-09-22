#!/bin/bash

for fullfile in /home/ubuntu/polyomica/nostaeva/Back_Pain_2021/validation_models/bp_adj/*.snpRes; do
echo ${fullfile}

dirname=$(dirname -- "$fullfile")
filename=$(basename -- "$fullfile")
extension="${filename##*.}"
filename="${filename%.*}"

echo gwas_id	chr	ea	ra	eaf	beta	se > ${dirname}/${filename}_genome.csv
awk '{print $2" "$3" "$5" "$6" "$7" "$8" "$9}' ${fullfile} | grep -v Name >> ${dirname}/${filename}_genome.csv
done
