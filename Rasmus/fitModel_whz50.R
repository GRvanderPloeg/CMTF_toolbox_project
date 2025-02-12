library(tidyverse)
library(ggplot2)
library(stringr)
library(parafac4microbiome)
library(CMTFtoolbox)

Y = Jakobsen2025$subjectMeta_WHZ$whz.6m
Ycnt = Y - mean(Y)
Ynorm = Ycnt / norm(Ycnt, "2")

model = CMTFtoolbox::acmtfr_opt(Jakobsen2025$Zwhz, Ynorm, numComponents=2, pi=0.50, nstart=100, numCores=parallel::detectCores(), method="L-BFGS", allOutput=TRUE)
saveRDS(model, "./ACMTFR_model_whz50.RDS")
