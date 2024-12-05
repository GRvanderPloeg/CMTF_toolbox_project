library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)

set.seed(123)

# Load data
X1_final = readRDS("./20241203_X1_Y_inside.RDS")
X2_final = readRDS("./20241203_X2_Y_inside.RDS")
X3_final = readRDS("./20241203_X3_Y_inside.RDS")
Y_final = readRDS("./20241203_Y_Y_inside.RDS")

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=TRUE)

# Prepare sim settings
betas = 1 * 10^seq(0, -10, length.out=100)
pis   = seq(0, 1, length.out=10)
simSettings = cbind(rep(betas, each=10), rep(pis, 10))

# Run model
cl = parallel::makeCluster(32)
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:length(betas)) %dopar% {
  betaValue = simSettings[i,1]
  piValue = simSettings[i,2]
  model=CMTFtoolbox::acmtf_optr(Z,as.matrix(Y_final@data),numComponents=7,initialization="nvec",beta=rep(betasValue,3),pi=piValue,abs_tol=1e-10,rel_tol=1e-10,nstart=1)
}
parallel::stopCluster(cl)

# Save model
saveRDS(models, "20241203_ACMTFR_CV_Y_inside.RDS")
saveRDS(simSettings, "20241203_ACMTFR_CV_params_Y_inside.RDS")
