#!/bin/bash
#
#SBATCH --job-name=sim1b_acmtf_cv
#SBATCH --output=res_Sim1b_ACMTF_CV.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=1-00:30:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=g.r.ploeg@uva.nl

date
Rscript Sim1b_ACMTF_CV_Y_inside.R
date
