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
datasets = list(X1_final, X2_final, X3_final, matrix(Y_final@data))
modes = list(c(1,2,3), c(1,4,5), c(1,6,7), c(1,8))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=TRUE)

# Prepare sim settings
betas = c(1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1)
betas = rep(betas, each=100)

# Run model
cl = parallel::makeCluster(64)
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:length(betas)) %dopar% {
  model=CMTFtoolbox::acmtf_opt(Z,numComponents=7,initialization="random",beta=rep(betas[i],4),abs_tol=1e-10,rel_tol=1e-10,nstart=1)
}
parallel::stopCluster(cl)

# Save model
saveRDS(models, "Sim2_ACMTF_YAS_models_Y_inside.RDS")
