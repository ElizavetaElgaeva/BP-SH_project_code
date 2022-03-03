# Aim of this script is to calculate phenotypes for SGIT,
# CBP UGIT and unique genetic impact of CBP trait
# for all Europeans from UKBB cohort

library(data.table)

setwd("/mnt/polyomica/projects/bp-sh/")

# Read a table with all pain phenotypes
dat <- fread("/home/freydin/Pheno/mv_mar18_chron.phen.txt", data.table = F)
dim(dat)
colnames(dat)
table(dat$ethn5)
#dat <- dat[!is.na(dat$ethn5) & dat$ethn5 == "White", ]
#dim(dat)
table(is.na(dat$back))
dat <- dat[!is.na(dat$back), ]
table(is.na(dat$neck)) # all FALSE
table(is.na(dat$knee)) # all FALSE
table(is.na(dat$hip)) # all FALSE
table(is.na(dat$head)) # all FALSE
table(is.na(dat$stom)) # all FALSE
dim(dat)

# Linear combination coefficents for SGIT
traits <- c("hip", "back", "neck", "knee", "head", "stom") # trait order for linear combination 

sgit <- read.table('./data/01_sh/alphas.txt', row.names = 1)
sgit <- as.numeric(sgit[2, ])

# Linear combination coefficents for CBP UGIT
n_traits <- c(1:length(traits))

source("./src/shared_heredity/00_core_functions/gcov_for_linear_comb_with_i_trait.R")
source("./src/shared_heredity/00_core_functions/heritability_of_linear_combination.R")

phem <- read.table('./data/01_sh/pheno_corr_matrix.txt', row.names = 1,  check.names = F)
gcov <- read.table('./data/01_sh/gene_cov_matrix.txt', row.names = 1,  check.names = F)

slope <- sapply(n_traits, function(x) cov_gi_alpha(a = sgit, i = x, covm = as.matrix(gcov))/H2(sgit, covm = gcov, phem = phem))
position <- diag(length(traits))
ugit <- lapply(n_traits, function(x) return(position[x, ] - sgit*slope[x]))
ugit <- ugit[[2]]

# Linear combination coefficents for unique genetic impact of CBP trait
unigit <- c(-0.071018123, 0.854891273, -0.790903139, 0.017435276, 0.001350986, -0.067428072) # "Hip", "Back", "Neck", "Knee", "Head", "Stomach"

# Standardize phenotypes
std_hip <- dat$hip/sd(dat$hip)
std_back <- dat$back/sd(dat$back)
std_neck <- dat$neck/sd(dat$neck)
std_knee <- dat$knee/sd(dat$knee)
std_head <- dat$head/sd(dat$head)
std_stom <- dat$stom/sd(dat$stom)

# Compute new phenotypes
dat$sgit <- std_hip*sgit[1] + std_back*sgit[2] + std_neck*sgit[3] + std_knee*sgit[4] + std_head*sgit[5] + std_stom*sgit[6]
dat$cbp_ugit <- std_hip*ugit[1] + std_back*ugit[2] + std_neck*ugit[3] + std_knee*ugit[4] + std_head*ugit[5] + std_stom*ugit[6]
dat$cbp_unique <- std_hip*unigit[1] + std_back*unigit[2] + std_neck*unigit[3] + std_knee*unigit[4] + std_head*unigit[5] + std_stom*unigit[6]

var(dat$sgit)
#[1] 1.018469
var(dat$cbp_ugit)
#[1] 0.5191914
var(dat$cbp_unique)
#[1] 0.9780631
cor(dat$sgit, dat$cbp_ugit)
#[1] -0.01658693

fwrite(dat, "./data/12_prs/phenotypic_data.txt", dec = ".", sep = "\t")

