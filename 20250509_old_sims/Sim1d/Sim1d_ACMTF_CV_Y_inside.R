library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(parallel)
library(doParallel)
library(foreach)

set.seed(123)

numSubjects = c(50, 200, 500)
numFeatures = c(100, 1000)
numTimepoints = c(10, 50, 100)
betas = c(1e-5, 1e-4, 1e-3, 1e-2)
l = list("numSubjects"=numSubjects, "numFeatures"=numFeatures, "numTimepoints"=numTimepoints, "beta"=betas)
simSettings = do.call(expand.grid, l)
simSettings = do.call(rbind, replicate(100, simSettings, simplify=FALSE))

# Run model
cl = parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:nrow(simSettings)) %dopar% {

  numSubjects = simSettings[i,1]
  numFeatures = simSettings[i,2]
  numTimepoints = simSettings[i,3]
  beta = simSettings[i,4]

  # Load data
  X1_final = readRDS(paste0("./Sim1d_X1_Y_inside", "_", numSubjects, "_", numFeatures, "_", numTimepoints, ".RDS"))
  X2_final = readRDS(paste0("./Sim1d_X2_Y_inside", "_", numSubjects, "_", numFeatures, "_", numTimepoints, ".RDS"))
  X3_final = readRDS(paste0("./Sim1d_X3_Y_inside", "_", numSubjects, "_", numFeatures, "_", numTimepoints, ".RDS"))

  # Process data
  datasets = list(X1_final, X2_final, X3_final)
  modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
  Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

  model=CMTFtoolbox::acmtf_opt(Z,numComponents=7,initialization="random",beta=rep(beta,3),abs_tol=1e-10,rel_tol=1e-10,nstart=1)
}
parallel::stopCluster(cl)

# Save model and parameters
saveRDS(models, paste0("Sim1d_ACMTF_CV_models_Y_inside.RDS"))
saveRDS(simSettings, paste0("Sim1d_ACMTF_CV_params_Y_inside.RDS"))
