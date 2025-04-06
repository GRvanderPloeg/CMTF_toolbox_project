library(tidyverse)
library(ggplot2)
library(stringr)
library(parafac4microbiome)
library(CMTFtoolbox)

result = ncrossreg(Jakobsen2025$Zwhz, Jakobsen2025$subjectMeta_WHZ$whz.6m, maxNumComponents=5, pi=0.25, nstart=10, numCores=parallel::detectCores(), method="L-BFGS", cvFolds=10)
saveRDS(result, paste0("./CV_manual_whz25.RDS"))



