This folder contains all code extracting shared heredity from 6 chronic pain traits in Discovery and Replication studies. Core functions used in this pipeline obtained from https://github.com/Sodbo/shared_heredity repository.

## 00a_start.sh
This script sources 00b_pipeline_for_matrix_calculation.sh script and sets the command line variables for further calculations (only original traits are mentioned). 

## 00b_pipeline_for_matrix_calculation.sh
This script incorporates a pipeline for calculation of correlation and covariance matrices for original traits and estimation of alpha coeffitients for linear combination (see steps 01 - 06). The script uses command line variables.

## 01_pheno_corr.sh
This script calculates phenotypic correlations between traits using data from Discovery study. The script is writted for GWAS-Map database and it uses command line variables.

## 02_convert_long_to_wide_form.R
This script converts phenotypic correlations matrix from long to wide format. The script uses command line variables.

## 03_gene_corr.sh
This script calculates genetic correlations between traits using data from Discovery study. The script is writted for GWAS-Map database and it uses command line variables.

## 04_gene_corr_to_matrices.R
This script calculates genetic correlations, covariances and se of genetic covariances between traits and builds matrices (dased on data from Discovery study). The script uses command line variables.

## 05_alpha_coefficients.R
This script calculates alpha coeffitients and weights for each trait in linear combination using data from Discovery study. The script uses command line variables.

## 05a_alpha_coefficients_five_traits.R
This script calculates alpha coeffitients in linear combination of five pain traits (Head pain was excluded) and estimates the expected genetic correlation between the shared heredity for five and six pain traits. The idea of these manipulations was to check for the possible bias due to Head pain inclusion in the analysis.  

## 06_alpha_CI_estimation.R
This script estimates confident intervals for alpha coeffitients.

## 07_linear_combination_sh.R
This script calculates summary statistics for shared heredity in Discovery study.

## 07a_validity_check_sh.R
This script provides analytical estimation of heritability of shared heredity and its genetic correlations with original traits using core functions (Discovery study).

## 07b_validity_check_sh.sh
This script estimates heritability of shared heredity and its genetic correlations with original traits using LD Score regression implemented in GWAS-Map database (Discovery study).

## vim 07c_validity_check_sh_gip1_rg.R
This script estimates expected genetic correlation between shared heredity of six pain traits from this study and GIP1 of four musculosceletal pain traits from our prevoius study (https://doi.org/10.1038/s42003-020-1051-9).

## 08_linear_combination_sh_aa_repl.R
This script calculates summary statistics for shared heredity in AA Replication study (replication cohort of African ancestry).

## 08a_validity_check_sh_aa_repl.sh
This script estimates heritability of shared heredity and its genetic correlations with original traits using LD Score regression implemented in GWAS-Map database (AA Replication study).

## 09_linear_combination_sh_ea_repl.R
This script calculates summary statistics for shared heredity in EA Replication study (replication cohort of European ancestry).

## 09a_validity_check_sh_ea_repl.sh
This script estimates heritability of shared heredity and its genetic correlations with original traits using LD Score regression implemented in GWAS-Map database (EA Replication study).

## 10_linear_combination_sh_sa_repl.R
This script calculates summary statistics for shared heredity in SA Replication study (replication cohort of South Asian ancestry).

## 10a_validity_check_sh_sa_repl.sh
This script estimates heritability of shared heredity and its genetic correlations with original traits using LD Score regression implemented in GWAS-Map database (SA Replication study).
