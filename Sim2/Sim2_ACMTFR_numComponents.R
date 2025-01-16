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
X1_final = readRDS("./Sim2_X1_Y_inside.RDS")
X2_final = readRDS("./Sim2_X2_Y_inside.RDS")
X3_final = readRDS("./Sim2_X3_Y_inside.RDS")
Y_final = readRDS("./Sim2_Y_Y_inside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

plot = CMTFtoolbox::ncrossreg(Z, as.matrix(Y_final@data), maxNumComponents=10, nstart=10, numCores=parallel::detectCores())
saveRDS(plot, "CV_numComponents.RDS")
