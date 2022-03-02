This folder contains scripts for analysis of genetic correlations between shared heredity (European meta-analysis) and other complex traits.

## 01_calculation_of_rg.sh
This script is for calculation of genetic correlations between SH_EUR_MA and all of the non eQTL GWASes with sample size>10000 (h2>0.001 and Z-score(h2)>3)

## 02_extract_rg_info
This script to extract correlations from GWAS-MAP to file.

## 03_select_ids_for_squared_rg_matrix.R
This script selects GWAS ids from GWAS-MAP database to create a squared matriz of genetic correlations.

## 04_
This scripts count upper-corner diagonal of gene corr matrix using LD Score. To run the script it is necessary to set a filename with gwas_ids (should be listed in column) as the command line argument.

## 05_
This scripts download the resulting matrix of genetic correlations from the GWAS-MAP.

## 06_corr_matrix_square_form.R
This script transforms the genetic correlation matrix to a square form.

## 07_hclustering.R
This script performs hierarchical clustering of the traits statistically significantly genetically correlated with either SH or BP-SH with |rg| > 0.25.

## 08_heatmap.R
This script creates a heatmap visualisation of the genetic correlation matrix based on clusters defined at the previous step.
