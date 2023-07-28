#!/bin/bash
#
#SBATCH --job-name=run
#SBATCH --output=comb.smp.log

srun R --vanilla < imp.sumFr_sgit.R
srun R --vanilla < imp.comb_meth.R
srun R --vanilla < imp.comb_sample.R
