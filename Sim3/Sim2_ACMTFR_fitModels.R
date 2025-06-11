library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(parallel)
library(doParallel)
library(foreach)

set.seed(123)

# Settings
gammas = c(0.001, 0.01, 0.1, 0.25, 0.50, 0.75, 1)

for(i in 1:length(gammas)){
  gamma = gammas[i]

  # Load data
  X1_final = readRDS(paste0("./Sim1_X1_", gamma, ".RDS"))
  X2_final = readRDS(paste0("./Sim1_X2_", gamma, ".RDS"))
  X3_final = readRDS(paste0("./Sim1_X3_", gamma, ".RDS"))
  Y_final = as.matrix(readRDS("./Sim1_Y.RDS"))

  # Process data
  datasets = list(X1_final, X2_final, X3_final)
  modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
  Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

  # Run ACMTF-R across pi values
  result = acmtfr_opt(Z, Y_final, 7, pi=0.0, method="L-BFGS", nstart=10, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim1_ACMTFR_model_00_", gamma, ".RDS"))

  result = acmtfr_opt(Z, Y_final, 7, pi=0.1, method="L-BFGS", nstart=10, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim1_ACMTFR_model_01_", gamma, ".RDS"))

  result = acmtfr_opt(Z, Y_final, 7, pi=0.2, method="L-BFGS", nstart=10, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim1_ACMTFR_model_02_", gamma, ".RDS"))

  result = acmtfr_opt(Z, Y_final, 7, pi=0.3, method="L-BFGS", nstart=10, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim1_ACMTFR_model_03_", gamma, ".RDS"))

  result = acmtfr_opt(Z, Y_final, 7, pi=0.4, method="L-BFGS", nstart=10, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim1_ACMTFR_model_04_", gamma, ".RDS"))

  result = acmtfr_opt(Z, Y_final, 7, pi=0.5, method="L-BFGS", nstart=10, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim1_ACMTFR_model_05_", gamma, ".RDS"))

  result = acmtfr_opt(Z, Y_final, 7, pi=0.6, method="L-BFGS", nstart=10, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim1_ACMTFR_model_06_", gamma, ".RDS"))

  result = acmtfr_opt(Z, Y_final, 7, pi=0.7, method="L-BFGS", nstart=10, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim1_ACMTFR_model_07_", gamma, ".RDS"))

  result = acmtfr_opt(Z, Y_final, 7, pi=0.8, method="L-BFGS", nstart=10, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim1_ACMTFR_model_08_", gamma, ".RDS"))

  result = acmtfr_opt(Z, Y_final, 7, pi=0.9, method="L-BFGS", nstart=10, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim1_ACMTFR_model_09_", gamma, ".RDS"))

  result = acmtfr_opt(Z, Y_final, 7, pi=1.0, method="L-BFGS", nstart=10, numCores=parallel::detectCores(), allOutput=TRUE)
  saveRDS(result, paste0("./Sim1_ACMTFR_model_10_", gamma, ".RDS"))

}
