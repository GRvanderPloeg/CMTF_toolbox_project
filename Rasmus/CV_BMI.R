library(tidyverse)
library(ggplot2)
library(stringr)
library(parafac4microbiome)
library(CMTFtoolbox)

pi_values = c(0.5, 0.75, 0.9, 0.95, 1)
for(i in 1:length(pi_values)){
  pi = pi_values[i]
  print(pi)

  Y = Jakobsen2025$homogenizedSubjectMetadata$BMI
  Ycnt = Y - mean(Y)

  result = ncrossreg(Jakobsen2025$Z, Ycnt, maxNumComponents=10, pi=pi, nstart=50, numCores=parallel::detectCores(), method="L-BFGS", cvFolds=10)
  saveRDS(result, paste0("./CV_BMI_", pi, ".RDS"))
}


