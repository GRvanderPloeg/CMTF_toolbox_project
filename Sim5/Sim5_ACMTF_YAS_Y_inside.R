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
datasets = list(X1_final, X2_final, X3_final, matrix(Y_final@data))
modes = list(c(1,2,3), c(1,4,5), c(1,6,7), c(1,8))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=TRUE)

# Prepare sim settings
betas = c(1e-3)
betas = rep(betas, each=100)

# Run model
cl = parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:length(betas)) %dopar% {
  model=CMTFtoolbox::acmtf_opt(Z,numComponents=7,initialization="random",beta=rep(betas[i],4),abs_tol=1e-10,rel_tol=1e-10,nstart=1)
}
parallel::stopCluster(cl)

# Save model
saveRDS(models, "Sim5_ACMTF_YAS_models_Y_inside.RDS")
saveRDS(betas, "Sim5_ACMTF_YAS_params_Y_inside.RDS")
