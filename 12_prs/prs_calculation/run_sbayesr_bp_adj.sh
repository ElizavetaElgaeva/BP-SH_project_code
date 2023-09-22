#!/bin/bash

path_to_data="/home/ubuntu/polyomica/nostaeva/Back_Pain_2021/data/BP_adj_disc_output_done.ma"
path_to_out="/home/ubuntu/polyomica/nostaeva/Back_Pain_2021/validation_models/bp_adj"
path_to_ld="/home/ubuntu/polyomica/nostaeva/Back_Pain_2021/process_data/list_ld_matrix.txt"
version="gctb_2.03beta_Linux"
flag="BP_adj"

/home/ubuntu/code_folder/tools/${version}/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 5000 \
--burn-in 1000 \
--exclude-mhc \
--no-mcmc-bin \
--robust \
--out ${path_to_out}/${flag}_model_1 >> ${flag}_model_1.log

/home/ubuntu/code_folder/tools/${version}/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 5000 \
--burn-in 1000 \
--pi 0.7,0.3 \
--gamma 0,1 \
--exclude-mhc \
--no-mcmc-bin \
--robust \
--out ${path_to_out}/${flag}_model_2 >> ${flag}_model_2.log


/home/ubuntu/code_folder/tools/${version}/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 5000 \
--burn-in 1000 \
--pi 0.9,0.1 \
--gamma 0,1 \
--exclude-mhc \
--no-mcmc-bin \
--robust \
--out ${path_to_out}/${flag}_model_3 >> ${flag}_model_3.log

/home/ubuntu/code_folder/tools/${version}/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 5000 \
--burn-in 1000 \
--pi 0.97,0.03 \
--gamma 0,1 \
--exclude-mhc \
--no-mcmc-bin \
--robust \
--out ${path_to_out}/${flag}_model_4 >> ${flag}_model_4.log

/home/ubuntu/code_folder/tools/${version}/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 5000 \
--burn-in 1000 \
--pi 0.99,0.01 \
--gamma 0,1 \
--exclude-mhc \
--no-mcmc-bin \
--robust \
--out ${path_to_out}/${flag}_model_5 >> ${flag}_model_5.log

/home/ubuntu/code_folder/tools/${version}/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 5000 \
--burn-in 1000 \
--pi 0.997,0.003 \
--gamma 0,1 \
--exclude-mhc \
--no-mcmc-bin \
--robust \
--out ${path_to_out}/${flag}_model_6 >> ${flag}_model_6.log
