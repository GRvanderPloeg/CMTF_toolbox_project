library(tidyverse)
library(ggplot2)
library(stringr)
library(parafac4microbiome)
library(CMTFtoolbox)

pi_values = c(0.25, 0.5, 0.75, 0.9, 0.95, 1)
for(i in 1:length(pi_values)){
  pi = pi_values[i]
  print(pi)

  result = ncrossreg(Jakobsen2025$Zbmi, Jakobsen2025$subjectMeta_BMI$BMI, maxNumComponents=5, pi=pi, nstart=10, numCores=parallel::detectCores(), method="L-BFGS", cvFolds=10)
  print(result)
  saveRDS(result, paste0("./CV_small_BMI_", pi, ".RDS"))
}


