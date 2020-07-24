This folder contains scripts for subtraction of Shared heredity trait from Back pain in Discovery and Replication cohorts. Core functions used in this pipeline were obtained from https://github.com/Sodbo/shared_heredity repository.

## 01_linear_combination.R
This script substracts shared heredity from back pain in Discovery study.

## 01a_validity_check.R
This script provides analytical estimation of heritability of traits without shared heredity and their genetic correlations with shared heredity using core functions.

## 01b_validity_check.sh
This script estimates heritability of back pain without shared heredity and its genetic correlation with shared heredity using LD Score regression implemented in GWAS-Map (Discovery study).

## 02_linear_combination_aa_repl.R
This script substracts shared heredity from back pain in AA Replication study (cohort of African ancestry).

## 03_linear_combination_ea_repl.R
This script substracts shared heredity from back pain in EA Replication study (cohort of European ancestry).

## 03a_validity_check_ea_repl.sh
This script estimates heritability of back pain without shared heredity and its genetic correlation with shared heredity using LD Score regression implemented in GWAS-Map (EA Replication study).

## 04_linear_combination_sa_repl.R
This script substracts shared heredity from back pain in SA Replication study (cohort of South Asian ancestry).
