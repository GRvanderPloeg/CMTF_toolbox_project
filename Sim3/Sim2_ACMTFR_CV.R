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

  # Run CV across pi values
  result = ACMTFR_modelSelection(datasets, modes, Y_final, maxNumComponents=10, pi=0.00, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
  saveRDS(result, paste0("./Sim1_ACMTFR_CV_00_", gamma, ".RDS"))

  result = ACMTFR_modelSelection(datasets, modes, Y_final, maxNumComponents=10, pi=0.10, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
  saveRDS(result, paste0("./Sim1_ACMTFR_CV_01_", gamma, ".RDS"))

  result = ACMTFR_modelSelection(datasets, modes, Y_final, maxNumComponents=10, pi=0.20, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
  saveRDS(result, paste0("./Sim1_ACMTFR_CV_02_", gamma, ".RDS"))

  result = ACMTFR_modelSelection(datasets, modes, Y_final, maxNumComponents=10, pi=0.30, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
  saveRDS(result, paste0("./Sim1_ACMTFR_CV_03_", gamma, ".RDS"))

  result = ACMTFR_modelSelection(datasets, modes, Y_final, maxNumComponents=10, pi=0.40, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
  saveRDS(result, paste0("./Sim1_ACMTFR_CV_04_", gamma, ".RDS"))

  result = ACMTFR_modelSelection(datasets, modes, Y_final, maxNumComponents=10, pi=0.50, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
  saveRDS(result, paste0("./Sim1_ACMTFR_CV_05_", gamma, ".RDS"))

  result = ACMTFR_modelSelection(datasets, modes, Y_final, maxNumComponents=10, pi=0.60, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
  saveRDS(result, paste0("./Sim1_ACMTFR_CV_06_", gamma, ".RDS"))

  result = ACMTFR_modelSelection(datasets, modes, Y_final, maxNumComponents=10, pi=0.70, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
  saveRDS(result, paste0("./Sim1_ACMTFR_CV_07_", gamma, ".RDS"))

  result = ACMTFR_modelSelection(datasets, modes, Y_final, maxNumComponents=10, pi=0.80, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
  saveRDS(result, paste0("./Sim1_ACMTFR_CV_08_", gamma, ".RDS"))

  result = ACMTFR_modelSelection(datasets, modes, Y_final, maxNumComponents=10, pi=0.90, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
  saveRDS(result, paste0("./Sim1_ACMTFR_CV_09_", gamma, ".RDS"))

  result = ACMTFR_modelSelection(datasets, modes, Y_final, maxNumComponents=10, pi=1.00, method="L-BFGS", nstart=10, cvFolds=10, numCores=parallel::detectCores())
  saveRDS(result, paste0("./Sim1_ACMTFR_CV_10_", gamma, ".RDS"))
}



