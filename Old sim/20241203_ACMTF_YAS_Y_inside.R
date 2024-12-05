library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)

set.seed(123)

# Load data
X1_final = readRDS("./20241203_X1_Y_inside.RDS")
X2_final = readRDS("./20241203_X2_Y_inside.RDS")
X3_final = readRDS("./20241203_X3_Y_inside.RDS")
Y_final = readRDS("./20241203_Y_Y_inside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final, matrix(Y_final@data))
modes = list(c(1,2,3), c(1,4,5), c(1,6,7), c(1,8))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=TRUE)

# Run model
models = acmtf_opt(Z, numComponents=7, beta=c(1e-3,1e-3,1e-3,1e-3), abs_tol=1e-10, rel_tol=1e-10, nstart=100, numCores=32, allOutput=TRUE)

# Save model
saveRDS(models, "20241203_ACMTF_YAS_models_Y_inside.RDS")
