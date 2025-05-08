library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(parallel)
library(doParallel)
library(foreach)

set.seed(123)
piValue = 0.5

# Load data
X1_final = readRDS("./Sim2b_X1_Y_inside.RDS")
X2_final = readRDS("./Sim2b_X2_Y_inside.RDS")
X3_final = readRDS("./Sim2b_X3_Y_inside.RDS")
Y_final = readRDS("./Sim2b_Y_Y_inside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))

# Run model CV
result=CMTFtoolbox::ACMTF_modelSelection(datasets,modes,maxNumComponents=10,nstart=10,cvFolds=10,method="L-BFGS",numCores=parallel::detectCores())
saveRDS(result, paste0("Sim2b_ACMTF_model_selection_result.RDS"))
