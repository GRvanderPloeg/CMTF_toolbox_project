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
X1_final = readRDS("./Sim3_X1_Y_outside.RDS")
X2_final = readRDS("./Sim3_X2_Y_outside.RDS")
X3_final = readRDS("./Sim3_X3_Y_outside.RDS")
Y_final = readRDS("./Sim3_Y_Y_outside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

# Prepare sim settings
betas = c(1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1)
pis   = seq(0, 1, length.out=6)
simSettings = cbind(rep(betas, each=6), rep(pis, 9))
simSettings = do.call(rbind, replicate(100, simSettings, simplify=FALSE))

# Run model
cl = parallel::makeCluster(64)
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:nrow(simSettings)) %dopar% {
  betaValue = simSettings[i,1]
  piValue = simSettings[i,2]
  model=CMTFtoolbox::acmtfr_opt(Z,as.matrix(Y_final@data),numComponents=7,initialization="random",beta=rep(betaValue,3),pi=piValue,abs_tol=1e-10,rel_tol=1e-10,nstart=1)
}
parallel::stopCluster(cl)

# Save model
saveRDS(models, "Sim3_ACMTFR_CV_Y_outside.RDS")
saveRDS(simSettings, "Sim3_ACMTFR_CV_params_Y_outside.RDS")
