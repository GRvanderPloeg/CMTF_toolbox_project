library(tidyverse)
library(ggplot2)
library(stringr)
library(parafac4microbiome)
library(CMTFtoolbox)

pi_values = c(0.25, 0.5, 0.75, 0.9, 0.95, 1)
for(i in 1:length(pi_values)){
  pi = pi_values[i]
  print(pi)

  result = ncrossreg(Jakobsen2025$Zwhz, Jakobsen2025$subjectMeta_WHZ$whz.6m, maxNumComponents=10, pi=pi, nstart=50, numCores=parallel::detectCores(), method="L-BFGS", cvFolds=10, normY=1)
  print(result)
  saveRDS(result, paste0("./CV_WHZ_", pi, ".RDS"))
}


