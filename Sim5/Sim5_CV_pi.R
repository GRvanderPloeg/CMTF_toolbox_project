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
X1_final = readRDS("./Sim5_X1_Y_inside.RDS")
X2_final = readRDS("./Sim5_X2_Y_inside.RDS")
X3_final = readRDS("./Sim5_X3_Y_inside.RDS")
Y_final = readRDS("./Sim5_Y_Y_inside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

# Calculate
result = CMTFtoolbox::cv_degeneracy(Z, as.matrix(Y_final@data), numComponents = 1:10, pis=c(0.25, 0.50, 0.75, 1), nstart = 10, numCores=32)

# Save model
saveRDS(result, "Sim5_ACMTFR_CV_pi.RDS")
