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
X1_final = readRDS("./Sim1_X1_Y_inside.RDS")
X2_final = readRDS("./Sim1_X2_Y_inside.RDS")
X3_final = readRDS("./Sim1_X3_Y_inside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

# Prepare sim settings
betas = c(1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1)
betas = rep(betas, each=1000)

# Run model
cl = parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:length(betas)) %dopar% {
  model=CMTFtoolbox::acmtf_opt(Z,numComponents=7,initialization="random",beta=rep(betas[i],3),abs_tol=1e-10,rel_tol=1e-10,nstart=1,method="L-BFGS")
}
parallel::stopCluster(cl)

# Save model and parameters
saveRDS(models, "Sim1_ACMTF_CV_models_Y_inside.RDS")
saveRDS(betas, "Sim1_ACMTF_CV_params_Y_inside.RDS")
