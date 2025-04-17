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
embeddings = c(1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.25, 0.5, 0.75, 1)
pis = c(0.2, 0.4, 0.6, 0.8, 0.85, 0.9, 0.95, 1)

l = list("embedding"=embeddings, "pi"=pis)
simSettings = do.call(expand.grid, l)
simSettings = do.call(rbind, replicate(100, simSettings, simplify=FALSE))

# Run model
cl = parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:nrow(simSettings)) %dopar% {
  embedding = simSettings[i,1]
  piValue = simSettings[i,2]

  # Load data
  X1_final = readRDS(paste0("./Sim5b_X1_", embedding, "_Y_inside.RDS"))
  X2_final = readRDS(paste0("./Sim5b_X2_", embedding, "_Y_inside.RDS"))
  X3_final = readRDS(paste0("./Sim5b_X3_", embedding, "_Y_inside.RDS"))
  Y_final = readRDS("./Sim5b_Y_Y_inside.RDS")

  # Process data
  datasets = list(X1_final, X2_final, X3_final)
  modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
  Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=FALSE)

  model=CMTFtoolbox::acmtfr_opt(Z,as.matrix(Y_final@data),numComponents=7,initialization="random",beta=rep(1e-3,3),pi=piValue,abs_tol=1e-10,rel_tol=1e-10,nstart=1)
}
parallel::stopCluster(cl)

# Save model
saveRDS(models, "Sim5b_ACMTFR_CV_Y_inside.RDS")
saveRDS(simSettings, "Sim5b_ACMTFR_CV_params_Y_inside.RDS")
