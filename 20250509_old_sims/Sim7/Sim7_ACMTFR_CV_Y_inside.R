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
X1_final = readRDS("./Sim7_X1_Y_inside.RDS")
X2_final = readRDS("./Sim7_X2_Y_inside.RDS")
X3_final = readRDS("./Sim7_X3_Y_inside.RDS")
Y_final = readRDS("./Sim7_Y_Y_inside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

# Prepare sim settings
mus = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
simSettings = rep(mus, each=100)

# Run model
cl = parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:length(simSettings)) %dopar% {
  muValue = simSettings[i]
  model=CMTFtoolbox::acmtfr_opt(Z,as.matrix(Y_final@data),numComponents=7,initialization="random",beta=rep(1e-3,3),pi=0.95,mu=muValue,method="L-BFGS",nstart=1)
}
parallel::stopCluster(cl)

# Save model
saveRDS(models, "Sim7_ACMTFR_CV_Y_inside.RDS")
saveRDS(simSettings, "Sim7_ACMTFR_CV_params_Y_inside.RDS")
