#!/bin/bash
#
#SBATCH --job-name=sim2_acmtfr_fit
#SBATCH --output=res_Sim2_ACMTFR_fitModels.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=7-02:30:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=g.r.ploeg@uva.nl
#SBATCH --mail-type=ALL

date
Rscript Sim2_ACMTFR_fitModels.R
date
