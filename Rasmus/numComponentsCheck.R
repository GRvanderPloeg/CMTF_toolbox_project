library(tidyverse)
library(ggplot2)
library(stringr)
library(parafac4microbiome)
library(CMTFtoolbox)


pi_values = c(0.25, 0.5, 0.75, 0.9, 1)
for(i in 1:length(pi_values)){
  pi = pi_values[i]
  print(pi)
  FMS_result = ncrossreg(Jakobsen2025$Z, Jakobsen2025$homogenizedSubjectMetadata$BMI, maxNumComponents=10, pi=pi, nstart=20, numCores=parallel::detectCores())
  saveRDS(FMS_result, paste0("./FMS_result_", pi, ".RDS"))
}


