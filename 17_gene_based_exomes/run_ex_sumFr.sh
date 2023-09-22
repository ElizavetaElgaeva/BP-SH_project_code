#!/bin/bash
#
#SBATCH --job-name=run
#SBATCH --output=run.log

srun R --vanilla < ex.prep4prep.R
srun R --vanilla < ex.sumFr_col.R
srun R --vanilla < ex.comb_meth_col.R
