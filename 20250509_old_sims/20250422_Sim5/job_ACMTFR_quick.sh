#!/bin/bash
#
#SBATCH --job-name=sim5_quick
#SBATCH --output=res_Sim5_quick.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=g.r.ploeg@uva.nl
#SBATCH --mail-type=ALL

date
Rscript Sim5_ACMTFR_CV_Y_inside_quick.R
date
