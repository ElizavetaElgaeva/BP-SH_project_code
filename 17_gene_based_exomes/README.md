This folder contains scripts for Gene-Based Analysis of SGIT.

For this analysis, summary statistics from EWAS and the matrices of correlations between genotypes of all variants within a gene are used.
Four variant annotations (LoF, LoF+missense, LoF+protein coding, and all intragenic SNPs) were analyzed. 
Variants were annotated using the Ensembl Variant Effect Predictor (VEP).
Previously, the genotypes of ultra-rare variants with MAC ≤ 10 were collapsed to a single variant.

## ex.prep4prep.R
Aim of this script is to prepare input score files for sumFregat R-package. 

## ex.sumFr_col.R
This script run sumFREGAT. Three methods are applied: SKAT-O, PCA, and ACAT-V.
Three variant annotations ara analyzed: protein coding (cod), protein non-coding (ncod), and nonsynonymous (nsyn) SNPs. 

## ex.comb_meth_col.R
Aim of this script is to combine results from different methods.

## run_ex_sumFr.sh
This script runs three others.
