#!/bin/bash
#
#SBATCH --job-name=sim1d_acmtf_cv
#SBATCH --output=res_Sim1d_ACMTF_CV.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=7-00:30:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=g.r.ploeg@uva.nl
#SBATCH --mail-type=ALL

date
Rscript Sim1d_ACMTF_CV_Y_inside.R
date
