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
X1_final = readRDS("./Sim1e_X1_Y_inside.RDS")
X2_final = readRDS("./Sim1e_X2_Y_inside.RDS")
X3_final = readRDS("./Sim1e_X3_Y_inside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

# Prepare sim settings
alphas = c(0.1, 0.5, 1, 2, 5, 10)
alphas = rep(alphas, each=100)

# Run model
cl = parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:length(alphas)) %dopar% {
  model=CMTFtoolbox::acmtf_opt(Z,numComponents=7,initialization="random",alpha=alphas[i],beta=rep(1e-3,3),abs_tol=1e-10,rel_tol=1e-10,nstart=1)
}
parallel::stopCluster(cl)

# Save model and parameters
saveRDS(models, "Sim1e_ACMTF_CV_models_Y_inside.RDS")
saveRDS(alphas, "Sim1e_ACMTF_CV_params_Y_inside.RDS")
