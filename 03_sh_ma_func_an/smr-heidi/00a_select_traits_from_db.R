# Aim of this script is to select traits (eQTLs and complex traits) for SMR/HEIDI analysis

library(data.table) 

setwd('/mnt/polyomica/projects/bp-sh/data/03_sh_ma_func_an/smr-heidi/')

# Read a file with descriptors from GWAS-MAP (from the September, 2021)
desc <- fread('./20210901_DESCRIPTORS.csv', data.table=F)

dim(desc)
# 1743250      44
colnames(desc)

table(desc$molecular_domain)
# complex trait Complex trait Complex_trait       Disease
#eQTL

# First filter by access restrictions
table(desc$"access_restrictions")
desc <- desc[desc$"access_restrictions" == '{public}', ]

# Split by trait type
eqtl <- desc[desc$molecular_domain == 'eQTL', ]
compl_tr <- desc[desc$molecular_domain == 'complex trait' | desc$molecular_domain == 'Complex trait' | desc$molecular_domain == 'Complex_trait' | desc$molecular_domain == 'Disease' | desc$molecular_domain == 'glycomics' | desc$molecular_domain == 'lipids' | desc$molecular_domain == 'metabolites' | desc$molecular_domain == 'proteins', ]

# Filter compl_tr
to_exc <- which(compl_tr$"trait_type" == 'binary' & as.numeric(compl_tr$n_cases) < 2000)
compl_tr <- compl_tr[-to_exc, ]

table(compl_tr$collection)
to_exc <- which(compl_tr$collection == 'SomaLogic_2017') # exclude SomaLogic 2017, as we keep the latest release
compl_tr <- compl_tr[-to_exc, ]
# other traits might be filtered manually

# Filter eqtl
 t(table(eqtl$tissue))
eqtl <- eqtl[eqtl$tissue == 'Muscle_Skeletal' | eqtl$tissue == 'Nerve_Tibial' | eqtl$tissue == 'Peripheral blood' | eqtl$tissue == 'Pituitary' | eqtl$tissue == 'Platelets (Peripheral blood)' | eqtl$tissue == 'Whole_Blood' |  eqtl$tissue == 'CD14+ monocytes (Peripheral blood)' |  eqtl$tissue == 'CD15+ granulocytes (Peripheral blood)' | eqtl$tissue == 'CD19+ B lymphocytes (Peripheral blood)' | eqtl$tissue == 'CD4+ T lymphocytes (Peripheral blood)' | eqtl$tissue == 'CD8+ T lymphocytes (Peripheral blood)' |  eqtl$tissue == 'Brain_Substantia_nigra' | eqtl$tissue == 'Brain_Spinal_cord_cervical_c-1' | eqtl$tissue == 'Brain_Putamen_basal_ganglia' | eqtl$tissue == 'Brain_Nucleus_accumbens_basal_ganglia' | eqtl$tissue == 'Brain_Hypothalamus' | eqtl$tissue == 'Brain_Hippocampus' | eqtl$tissue == 'Brain_Frontal_Cortex_BA9' | eqtl$tissue == 'Brain_Cortex' | eqtl$tissue == 'Brain_Cerebellum' | eqtl$tissue == 'Brain_Cerebellar_Hemisphere' | eqtl$tissue == 'Brain_Caudate_basal_ganglia' | eqtl$tissue == 'Brain_Anterior_cingulate_cortex_BA24' | eqtl$tissue == 'Brain_Amygdala' | eqtl$tissue == 'Dorsal root ganglia' | eqtl$tissue == 'blood', ]  

# Write
fwrite(compl_tr, file = './complex_traits.txt', sep = '\t')
fwrite(as.data.table(compl_tr$gwas_id), file = './complex_traits_ids.txt', sep = '\t', row.names = FALSE, col.names = FALSE)
fwrite(eqtl, file = './eqtls.txt', sep = '\t')
fwrite(as.data.table(eqtl$gwas_id), file = './eqtls_ids.txt', sep = '\t', row.names = FALSE, col.names = FALSE)



