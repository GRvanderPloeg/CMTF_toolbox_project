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
X1_final = readRDS("./Sim2b_X1_Y_inside.RDS")
X2_final = readRDS("./Sim2b_X2_Y_inside.RDS")
X3_final = readRDS("./Sim2b_X3_Y_inside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

# Prepare sim settings
noises = seq(0, 1, length.out=11)
noises = rep(noises, each=100)

# Run model
cl = parallel::makeCluster(64)
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:length(noises)) %dopar% {
  noiseOnY = noises[i]
  Y_final = readRDS(paste0("./Sim2b_Y_Y_inside_", noiseOnY, ".RDS"))
  model=CMTFtoolbox::acmtfr_opt(Z,as.matrix(Y_final@data),numComponents=7,initialization="random",beta=rep(1e-3,3),pi=0.5,abs_tol=1e-10,rel_tol=1e-10,nstart=1)
}
parallel::stopCluster(cl)

# Save model
saveRDS(models, "Sim2b_ACMTFR_CV_Y_inside.RDS")
saveRDS(noises, "Sim2b_ACMTFR_CV_params_Y_inside.RDS")
