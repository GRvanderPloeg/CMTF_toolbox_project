#!/bin/bash
#
#SBATCH --job-name=sim2b_acmtfr_modelSelection
#SBATCH --output=res_Sim2b_ACMTFR_modelSelection.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:30:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=g.r.ploeg@uva.nl
#SBATCH --mail-type=ALL

date
Rscript Sim2b_ACMTF_model_selection.R
date
