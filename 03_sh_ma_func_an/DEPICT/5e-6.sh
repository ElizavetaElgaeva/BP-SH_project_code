#!/bin/bash

### The aim of this script is to run DEPICT for meta-analysis cohort for SH COJO results with p-value treshold 5e-06

### Config file parameters:

#analysis_path: /mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/depict/5e-06

#gwas_summary_statistics_file:  /mnt/polyomica/projects/bp-sh/data/07_cojo/01_sh/sh_ma_eur_all_chr_5e-6.jma.cojo

#plink_executable: /usr/local/bin/plink

#association_pvalue_cutoff:  1

/home/common/projects/depict_software/DEPICT_v194/src/python/depict.py /mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/depict/5e-06/sh_ma_eur_all_chr.cfg
