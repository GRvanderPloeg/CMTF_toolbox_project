library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(parallel)
library(doParallel)
library(foreach)

set.seed(123)

# Prepare sim settings
# noises = seq(0, 10, length.out=51)
# noises = rep(noises, each=100)
noises = readRDS("./Sim1b_noise_Y_inside.RDS")
noises = rep(noises, each=100)

# Run model
cl = parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:length(noises)) %dopar% {

  noiseOnX = noises[i]

  # Load data
  X1_final = readRDS(paste0("./Sim1b_X1_Y_inside_", noiseOnX, ".RDS"))
  X2_final = readRDS(paste0("./Sim1b_X2_Y_inside_", noiseOnX, ".RDS"))
  X3_final = readRDS(paste0("./Sim1b_X3_Y_inside_", noiseOnX, ".RDS"))

  # Process data
  datasets = list(X1_final, X2_final, X3_final)
  modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
  Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

  # Run model
  model=CMTFtoolbox::acmtf_opt(Z,numComponents=7,initialization="random",beta=rep(1e-3,3),abs_tol=1e-10,rel_tol=1e-10,nstart=1)
}
parallel::stopCluster(cl)

# Save model and parameters
saveRDS(models, "Sim1b_ACMTF_CV_models_Y_inside.RDS")
saveRDS("")
