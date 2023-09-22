This folder contains scripts for Gene-Based Analysis of two traits, SGIT and UGIT.

For this analysis, summary statistics from GWAS and the matrices of correlations between genotypes of all variants within a gene are used.


## imp.prep4prep_sh.R
Aim of this script is to prepare input score files for sumFregat R-package. 

## imp.sumFr_sh.R
This script run sumFREGAT. Three methods are applied: SKAT-O, PCA, and ACAT-V.
Three variant annotations ara analyzed: protein coding (cod), protein non-coding (ncod), and nonsynonymous (nsyn) SNPs. 

## imp.comb_meth.R
Aim of this script is to combine results from different methods.

## run_imp_sumFr.sh
This script runs three others.
