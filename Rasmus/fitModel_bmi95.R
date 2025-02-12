library(tidyverse)
library(ggplot2)
library(stringr)
library(parafac4microbiome)
library(CMTFtoolbox)

Y = Jakobsen2025$subjectMeta_BMI$BMI
Ycnt = Y - mean(Y)
Ynorm = Ycnt / norm(Ycnt, "2")

model = CMTFtoolbox::acmtfr_opt(Jakobsen2025$Zbmi, Ynorm, numComponents=5, pi=0.95, nstart=100, numCores=parallel::detectCores(), method="L-BFGS", allOutput=TRUE)
saveRDS(model, "./ACMTFR_model_bmi95.RDS")
