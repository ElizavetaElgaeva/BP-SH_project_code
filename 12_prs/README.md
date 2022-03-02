This folder contains code for PRS analyses of chronic back pain (CBP), shared heredity (aka SGIT) and back pain adjusted for shared heredity (aka CBP UGIT).

## 01a_filter_icd10.R
This script performs validity checks of ICD-10 codes and filters the data to keep only individuals with back pain information and remove ICD-10 codes with low prevalence.

## 01b_filter_opcs.R
This script performs validity checks of OPCS codes and filters the data to keep only individuals with back pain information and remove OPCS codes with low prevalence.

## 01c_select_non-relatives.R
This script filters data on ICD-10 and OPCS codes to keep only non-relative individuals from the test sampe for further PRS analysis.

## 02_v2_standardize_prs.R
This script performs standartization of CBP, SGIT and CBP UGIT PRS.

## 03a_v2_icd10_vs_prs.R
This script runs generalized linear model analyses using ICD-10 codes as depending variables and sex, age, batch, 1-10 first principal components and PRS (CBP, SGIT, CBP UGIT) as covariates.

## 03b_v2_opcs_vs_prs.R
This script runs generalized linear model analyses using OPCS codes as depending variables and sex, age, batch, 1-10 first principal components and PRS (CBP, SGIT, CBP UGIT) as covariates.

## 04a_v2_filter_bp_prs_glm_results.R
This script filters results of generalized linear model analyses of CBP PRS to keep the statistically significant results only.

## 04b_v2_filter_sh_prs_glm_results.R
This script filters results of generalized linear model analyses of SGIT PRS to keep the statistically significant results only. 

## 04c_v2_filter_bp-sh_prs_glm_results.R
This script filters results of generalized linear model analyses of CBP UGIT PRS to keep the statistically significant results only.

## 05_hclust.R
This script performs clustering of the ICD-10 and OPCS codes.

## 06_v2_glm_table.R
This script creates a table with the results of generalized linear model analyses of PRS, ICD-10 and OPCS data.

## 07_v2_hclust.R 
This script performs clustering of the ICD-10 and OPCS codes using non-relatieves only.

## 07_v2_heatmap.R
This script visualizes the results of generalized linear model analyses of PRS, ICD-10 and OPCS data based on clusters defined at the previous step.

## 08_auc_models.R
This script estimates the AUC values of generalized linear models including PRS data.
