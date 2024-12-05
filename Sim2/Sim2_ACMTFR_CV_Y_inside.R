library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)

set.seed(123)

# Load data
X1_final = readRDS("./Sim2_X1_Y_inside.RDS")
X2_final = readRDS("./Sim2_X2_Y_inside.RDS")
X3_final = readRDS("./Sim2_X3_Y_inside.RDS")
Y_final = readRDS("./Sim2_Y_Y_inside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

# Prepare sim settings
betas = 1 * 10^seq(0, -10, length.out=11)
pis   = seq(0, 1, length.out=11)
simSettings = cbind(rep(betas, each=11), rep(pis, 11))

# Run model
cl = parallel::makeCluster(64)
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:nrow(simSettings)) %dopar% {
  betaValue = simSettings[i,1]
  piValue = simSettings[i,2]
  model=CMTFtoolbox::acmtfr_opt(Z,as.matrix(Y_final@data),numComponents=7,initialization="random",beta=rep(betaValue,3),pi=piValue,abs_tol=1e-10,rel_tol=1e-10,nstart=100)
}
parallel::stopCluster(cl)

# Save model
saveRDS(models, "Sim2_ACMTFR_CV_Y_inside.RDS")
saveRDS(simSettings, "Sim2_ACMTFR_CV_params_Y_inside.RDS")
