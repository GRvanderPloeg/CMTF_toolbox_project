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
X1_final = readRDS("./20241203_X1_Y_inside.RDS")
X2_final = readRDS("./20241203_X2_Y_inside.RDS")
X3_final = readRDS("./20241203_X3_Y_inside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=TRUE)

# Prepare sim settings
betas = 1 * 10^seq(0, -10, length.out=100)

# Run model
cl = parallel::makeCluster(32)
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:length(betas)) %dopar% {
  model=CMTFtoolbox::acmtf_opt(Z,numComponents=7,initialization="nvec",beta=rep(betas[i],3),abs_tol=1e-10,rel_tol=1e-10,nstart=1)
}
parallel::stopCluster(cl)

# Save model and parameters
saveRDS(models, "20241203_ACMTF_CV_models_Y_inside.RDS")
saveRDS(betas, "20241203_ACMTF_CV_params_Y_inside.RDS")
