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
noises = noises / (1 - noises)
noises[11] = 999
pis = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0)
simSettings = cbind(rep(noises, each=13), rep(pis, 11))
simSettings = do.call(rbind, replicate(100, simSettings, simplify=FALSE))

# Run model
cl = parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)

models = foreach::foreach(i=1:nrow(simSettings)) %dopar% {
  noiseOnY = simSettings[i,1]
  pi_value = simSettings[i,2]
  Y_final = readRDS(paste0("./Sim2b_Y_Y_inside_", noiseOnY, ".RDS"))
  model=CMTFtoolbox::acmtfr_opt(Z,as.matrix(Y_final@data),numComponents=7,initialization="random",beta=rep(1e-3,3),pi=pi_value,abs_tol=1e-10,rel_tol=1e-10,nstart=1)
}
parallel::stopCluster(cl)

# Save model
saveRDS(models, "Sim2b_ACMTFR_CV_Y_inside.RDS")
saveRDS(noises, "Sim2b_ACMTFR_CV_params_Y_inside.RDS")
