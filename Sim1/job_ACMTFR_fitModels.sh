#!/bin/bash
#
#SBATCH --job-name=sim1_acmtfr_fit
#SBATCH --output=res_Sim1_ACMTFR_fitModels.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=02:30:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=g.r.ploeg@uva.nl
#SBATCH --mail-type=ALL

date
Rscript Sim1_ACMTFR_fitModels.R
date
