#!/bin/bash
#
#SBATCH --job-name=sim5_acmtf_yas
#SBATCH --output=res_Sim5_ACMTF_YAS.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=1-00:30:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=g.r.ploeg@uva.nl
#SBATCH --mail-type=ALL

date
Rscript Sim5_ACMTF_YAS_Y_inside.R
date
