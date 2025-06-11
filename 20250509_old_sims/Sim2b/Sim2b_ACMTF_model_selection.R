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
Y_final = readRDS("./Sim2b_Y_Y_inside_1.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))

# Run model CV
result=CMTFtoolbox::ACMTFR_modelSelection(datasets,modes,as.matrix(Y_final@data),pi=piValue,nstart=10,cvFolds=10, normalize=FALSE)
saveRDS(result, paste0("Sim2b_ACMTFR_model_selection_result", piValue, ".RDS"))
