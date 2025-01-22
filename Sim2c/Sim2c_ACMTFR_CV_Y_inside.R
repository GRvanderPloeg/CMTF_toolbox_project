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
X1_final = readRDS("./Sim2c_X1_Y_inside.RDS")
X2_final = readRDS("./Sim2c_X2_Y_inside.RDS")
X3_final = readRDS("./Sim2c_X3_Y_inside.RDS")
Y_final = readRDS("./Sim2c_Y_Y_inside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

# Prepare sim settings
betas = c(1e-4, 1e-3, 1e-2, 1e-1)
pis   = c(0.4, 0.6, 0.8, 0.9, 0.95, 1)
simSettings = cbind(rep(betas, each=6), rep(pis, 4))
simSettings = do.call(rbind, replicate(100, simSettings, simplify=FALSE))

# Run model
cl = parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:nrow(simSettings)) %dopar% {
  betaValue = simSettings[i,1]
  piValue = simSettings[i,2]
  model=CMTFtoolbox::acmtfr_opt(Z,as.matrix(Y_final@data),numComponents=7,initialization="random",beta=rep(betaValue,3),pi=piValue,abs_tol=1e-10,rel_tol=1e-10,nstart=1)
}
parallel::stopCluster(cl)

# Save model
saveRDS(models, "Sim2c_ACMTFR_CV_Y_inside.RDS")
saveRDS(simSettings, "Sim2c_ACMTFR_CV_params_Y_inside.RDS")
