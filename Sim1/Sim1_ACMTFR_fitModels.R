library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(parallel)
library(doParallel)
library(foreach)

set.seed(123)

# Load data
X1_final = readRDS("./Sim1_X1.RDS")
X2_final = readRDS("./Sim1_X2.RDS")
X3_final = readRDS("./Sim1_X3.RDS")
Y_final = as.matrix(readRDS("./Sim1_Y.RDS"))

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

# Number of models
numModels = 100

# Run ACMTF-R across pi values

# Lower values were not stable

# result = acmtfr_opt(Z, Y_final, 7, pi=0.6, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
# saveRDS(result, "./Sim1_ACMTFR_model_06.RDS")
#
# result = acmtfr_opt(Z, Y_final, 7, pi=0.7, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
# saveRDS(result, "./Sim1_ACMTFR_model_07.RDS")

result = acmtfr_opt(Z, Y_final, 4, pi=0.8, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
saveRDS(result, "./Sim1_ACMTFR_model_08.RDS")

# result = acmtfr_opt(Z, Y_final, 7, pi=0.9, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
# saveRDS(result, "./Sim1_ACMTFR_model_09.RDS")

result = acmtfr_opt(Z, Y_final, 6, pi=1.0, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
saveRDS(result, "./Sim1_ACMTFR_model_10.RDS")
