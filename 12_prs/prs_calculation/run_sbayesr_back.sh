#!/bin/bash

path_to_data="/home/ubuntu/polyomica/nostaeva/Back_Pain_2021/data/Back_output_done.ma"
path_to_out="/home/ubuntu/polyomica/nostaeva/Back_Pain_2021/validation_models"
path_to_ld="/home/ubuntu/polyomica/nostaeva/Back_Pain_2021/process_data/list_ld_matrix.txt"

cd /home/ubuntu/polyomica/nostaeva/Back_Pain_2021
mkdir validation_models

/home/ubuntu/code_folder/tools/gctb2.02/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 5000 \
--burn-in 1000 \
--exclude-mhc \
--no-mcmc-bin \
--out ${path_to_out}/Back_5000

/home/ubuntu/code_folder/tools/gctb2.02/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 3000 \
--burn-in 1000 \
--exclude-mhc \
--no-mcmc-bin \
--out ${path_to_out}/Back_3000


/home/ubuntu/code_folder/tools/gctb2.02/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 5000 \
--burn-in 1000 \
--pi 0.9,0.1 \
--gamma 0,1 \
--exclude-mhc \
--no-mcmc-bin \
--out ${path_to_out}/Back_pi_0.9_0.1_5000

/home/ubuntu/code_folder/tools/gctb2.02/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 5000 \
--burn-in 1000 \
--pi 0.7,0.3 \
--gamma 0,1 \
--exclude-mhc \
--no-mcmc-bin \
--out ${path_to_out}/Back_pi_0.7_0.3_5000

/home/ubuntu/code_folder/tools/gctb2.02/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 5000 \
--burn-in 1000 \
--pi 0.5,0.5 \
--gamma 0,1 \
--exclude-mhc \
--no-mcmc-bin \
--out ${path_to_out}/Back_pi_0.5_0.5_5000

/home/ubuntu/code_folder/tools/gctb2.02/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 5000 \
--burn-in 1000 \
--pi 0.3,0.7 \
--gamma 0,1 \
--exclude-mhc \
--no-mcmc-bin \
--out ${path_to_out}/Back_pi_0.3_0.7_5000

/home/ubuntu/code_folder/tools/gctb2.02/gctb --sbayes R \
--mldm ${path_to_ld} \
--gwas-summary ${path_to_data} \
--chain-length 5000 \
--burn-in 1000 \
--pi 0.1,0.9 \
--gamma 0,1 \
--exclude-mhc \
--no-mcmc-bin \
--out ${path_to_out}/Back_pi_0.1_0.9_5000
