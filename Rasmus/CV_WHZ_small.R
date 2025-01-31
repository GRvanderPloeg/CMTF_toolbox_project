library(tidyverse)
library(ggplot2)
library(stringr)
library(parafac4microbiome)
library(CMTFtoolbox)

pi_values = c(0.5, 0.75, 0.9, 0.95, 1)
for(i in 1:length(pi_values)){
  pi = pi_values[i]
  print(pi)

  # Re-define Y and Z due to missing WHZ data
  Y = Jakobsen2025$homogenizedSubjectMetadata$whz.6m
  mask = !is.na(Y)
  Y = Y[mask]
  Ycnt = Y - mean(Y)

  datasets = lapply(Jakobsen2025$Z$object, FUN=function(x){x@data[mask,,]})
  modes = Jakobsen2025$Z$modes
  Z = setupCMTFdata(datasets, modes)

  result = ncrossreg(Z, Ycnt, maxNumComponents=5, pi=pi, nstart=10, numCores=parallel::detectCores(), method="L-BFGS")
  saveRDS(result, paste0("./CV_small_WHZ_", pi, ".RDS"))
}


