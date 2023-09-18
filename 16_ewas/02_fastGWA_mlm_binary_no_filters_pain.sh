# Aim of this script is to run EWAS on pain phenotype using exome genotypes before and after collapsing
# White Europeans, UKBB

########################################################## Exomes ##################################################################

# CKP

/home/common/projects/CAD_exomes/software/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
	--grm-sparse /home/common/DataStorage/UKBB/Project_59345/GRM/ukb_ea_sp \
	--fastGWA-mlm-binary \
	--pfile /home/common/DataStorage/UKBB/Project_59345/Exome/ukb23155/white_poly_qc/ukb23155_w_p_qc \
	--autosome \
	--threads 6 \
	--maf 0.000005 \
	--geno 0.02 \
	--pheno /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CKP_self-rep_eur_gcta.txt \
	--covar /home/common/projects/varicose_project/2023/pr59345_bcovar_sex_batch.txt \
	--qcovar /home/common/projects/varicose_project/2023/pr59345_qcovar_pc_age.txt \
	--out /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CKP/CKP_self-rep_eur_exomes

# CBP

/home/common/projects/CAD_exomes/software/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
        --grm-sparse /home/common/DataStorage/UKBB/Project_59345/GRM/ukb_ea_sp \
        --fastGWA-mlm-binary \
        --pfile /home/common/DataStorage/UKBB/Project_59345/Exome/ukb23155/white_poly_qc/ukb23155_w_p_qc \
        --autosome \
        --threads 6 \
        --maf 0.000005 \
        --geno 0.02 \
	--pheno /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CBP_self-rep_eur_gcta_v2.txt \
	--covar /home/common/projects/varicose_project/2023/pr59345_bcovar_sex_batch.txt \
	--qcovar /home/common/projects/varicose_project/2023/pr59345_qcovar_pc_age.txt \
	--out /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CBP/CBP_self-rep_eur_exomes_v2

# CNP

/home/common/projects/CAD_exomes/software/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
        --grm-sparse /home/common/DataStorage/UKBB/Project_59345/GRM/ukb_ea_sp \
        --fastGWA-mlm-binary \
        --pfile /home/common/DataStorage/UKBB/Project_59345/Exome/ukb23155/white_poly_qc/ukb23155_w_p_qc \
        --autosome \
        --threads 6 \
        --maf 0.000005 \
        --geno 0.02 \
        --pheno /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CNP_self-rep_eur_gcta.txt \
        --covar /home/common/projects/varicose_project/2023/pr59345_bcovar_sex_batch.txt \
        --qcovar /home/common/projects/varicose_project/2023/pr59345_qcovar_pc_age.txt \
	--out /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CNP/CNP_self-rep_eur_exomes


# CHP

/home/common/projects/CAD_exomes/software/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
        --grm-sparse /home/common/DataStorage/UKBB/Project_59345/GRM/ukb_ea_sp \
        --fastGWA-mlm-binary \
        --pfile /home/common/DataStorage/UKBB/Project_59345/Exome/ukb23155/white_poly_qc/ukb23155_w_p_qc \
        --autosome \
	--threads 6 \
	--maf 0.000005 \
	--geno 0.02 \
	--pheno /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CHP_self-rep_eur_gcta.txt \
        --covar /home/common/projects/varicose_project/2023/pr59345_bcovar_sex_batch.txt \
        --qcovar /home/common/projects/varicose_project/2023/pr59345_qcovar_pc_age.txt \
        --out /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CHP/CHP_self-rep_eur_exomes



################################################## Collapsing ###########################################################################

cd /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS

for set in 'set1.mac10' 'set1.maf0.01' 'set2.mac10' 'set3.mac10' 'set4.mac10'
do

# CKP

	/home/common/projects/CAD_exomes/software/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
		--grm-sparse /home/common/DataStorage/UKBB/Project_59345/GRM/ukb_ea_sp \
		--fastGWA-mlm-binary \
		--bfile /home/common/DataStorage/UKBB/Project_59345/Exome/collapsing/bbf_r18219/all_chr_${set} \
		--autosome \
		--threads 6 \
		--maf 0.000005 \
		--geno 0.02 \
		--pheno /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CKP_self-rep_eur_gcta.txt \
		--covar /home/common/projects/varicose_project/2023/pr59345_bcovar_sex_batch.txt \
		--qcovar /home/common/projects/varicose_project/2023/pr59345_qcovar_pc_age.txt \
		--out /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CKP/CKP_self-rep_eur_exomes_${set}_collaps


# CBP

	/home/common/projects/CAD_exomes/software/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
		--grm-sparse /home/common/DataStorage/UKBB/Project_59345/GRM/ukb_ea_sp \
		--fastGWA-mlm-binary \
		--bfile /home/common/DataStorage/UKBB/Project_59345/Exome/collapsing/bbf_r18219/all_chr_${set} \
		--autosome \
		--threads 6 \
		--maf 0.000005 \
		--geno 0.02 \
		--pheno /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CBP_self-rep_eur_gcta_v2.txt \
		--covar /home/common/projects/varicose_project/2023/pr59345_bcovar_sex_batch.txt \
		--qcovar /home/common/projects/varicose_project/2023/pr59345_qcovar_pc_age.txt \
		--out /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CBP/CBP_self-rep_eur_exomes_${set}_collaps_v2


# CNP

	/home/common/projects/CAD_exomes/software/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
		--grm-sparse /home/common/DataStorage/UKBB/Project_59345/GRM/ukb_ea_sp \
		--fastGWA-mlm-binary \
		--bfile /home/common/DataStorage/UKBB/Project_59345/Exome/collapsing/bbf_r18219/all_chr_${set} \
		--autosome \
		--threads 6 \
		--maf 0.000005 \
		--geno 0.02 \
		--pheno /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CNP_self-rep_eur_gcta.txt \
		--covar /home/common/projects/varicose_project/2023/pr59345_bcovar_sex_batch.txt \
		--qcovar /home/common/projects/varicose_project/2023/pr59345_qcovar_pc_age.txt \
		--out /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CNP/CNP_self-rep_eur_exomes_${set}_collaps


# CHP

	/home/common/projects/CAD_exomes/software/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
	        --grm-sparse /home/common/DataStorage/UKBB/Project_59345/GRM/ukb_ea_sp \
		--fastGWA-mlm-binary \
		--bfile /home/common/DataStorage/UKBB/Project_59345/Exome/collapsing/bbf_r18219/all_chr_${set} \
		--autosome \
		--threads 6 \
		--maf 0.000005 \
		--geno 0.02 \
		--pheno /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/CHP_self-rep_eur_gcta.txt \
		--covar /home/common/projects/varicose_project/2023/pr59345_bcovar_sex_batch.txt \
		--qcovar /home/common/projects/varicose_project/2023/pr59345_qcovar_pc_age.txt \
		--out /home/common/projects/pain_project/fastGWA_mlm_binary_GWAS/no_filters/exomes/CHP/CHP_self-rep_eur_exomes_${set}_collaps



