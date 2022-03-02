This folder contains code generating and visualizing correlation matrices between original pain traits, shared heredity and back pain adjusted for shared heredity.

## 00a_start.sh
This code passes the gwas_id of the traits of interest and runs script 00b.

## 00b_pipeline_for_matrix_calculation.sh
This code runs scripts 01-04 using arguments passed in 00a.

## 01_pheno_corr.sh
This code calculates phenotypic correlation matrix using GWAS-Map.

## 02_convert_long_to_wide_form.R
This code transforms phenotypic correlation matrix to wide format.

## 03_gene_corr.sh
This code calculates pairwise genotypic correlations using GWAS-Map.

## 04_gene_corr_to_matrices.R
This code creates genotypic correlation matrix.

## 05_matrix_visualisation.R
This code visualizes phenotypic and genotypic correlation matrices as a heatmap.
