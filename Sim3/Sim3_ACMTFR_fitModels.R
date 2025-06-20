library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(parallel)
library(doParallel)
library(foreach)

set.seed(456)

# Settings
noises = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999)
numModels = 100
numComponents = 8

# Load Y once
Y_final = as.matrix(readRDS("./Sim3_Y.RDS"))

for(i in 1:length(noises)){

  noiseOnX = noises[i]

  # Load data
  X1_final = readRDS(paste0("./Sim3_X1_", noiseOnX, ".RDS"))
  X2_final = readRDS(paste0("./Sim3_X2_", noiseOnX, ".RDS"))
  X3_final = readRDS(paste0("./Sim3_X3_", noiseOnX, ".RDS"))

  # Process data
  datasets = list(X1_final, X2_final, X3_final)
  modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
  Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=TRUE)

  # Run ACMTF-R across pi values
  result = acmtfr_opt(Z, Y_final, numComponents, pi=0.1, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim3_ACMTFR_model_01_", noiseOnX, ".RDS"))

  result = acmtfr_opt(Z, Y_final, numComponents, pi=0.2, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim3_ACMTFR_model_02_", noiseOnX, ".RDS"))

  result = acmtfr_opt(Z, Y_final, numComponents, pi=0.3, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim3_ACMTFR_model_03_", noiseOnX, ".RDS"))

  result = acmtfr_opt(Z, Y_final, numComponents, pi=0.4, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim3_ACMTFR_model_04_", noiseOnX, ".RDS"))

  result = acmtfr_opt(Z, Y_final, numComponents, pi=0.5, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim3_ACMTFR_model_05_", noiseOnX, ".RDS"))

  result = acmtfr_opt(Z, Y_final, numComponents, pi=0.6, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim3_ACMTFR_model_06_", noiseOnX, ".RDS"))

  result = acmtfr_opt(Z, Y_final, numComponents, pi=0.7, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim3_ACMTFR_model_07_", noiseOnX, ".RDS"))

  result = acmtfr_opt(Z, Y_final, numComponents, pi=0.8, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim3_ACMTFR_model_08_", noiseOnX, ".RDS"))

  result = acmtfr_opt(Z, Y_final, numComponents, pi=0.9, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim3_ACMTFR_model_09_", noiseOnX, ".RDS"))

  result = acmtfr_opt(Z, Y_final, numComponents, pi=1.0, method="L-BFGS", nstart=numModels, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim3_ACMTFR_model_10_", noiseOnX, ".RDS"))

}


