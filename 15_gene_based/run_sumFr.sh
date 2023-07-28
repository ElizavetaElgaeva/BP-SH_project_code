#!/bin/bash
#
#SBATCH --job-name=run
#SBATCH --output=gba.log
srun R --vanilla < imp.prep4prep_sh.R
srun R --vanilla < imp.sumFr_sh.R
srun R --vanilla < imp.comb_meth.R
