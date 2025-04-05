library(tidyverse)
library(parafac4microbiome)
library(NPLStoolbox)
library(CMTFtoolbox)

homogenizedSamples = intersect(intersect(NPLStoolbox::Jakobsen2025$faeces$mode1$subject, NPLStoolbox::Jakobsen2025$milkMicrobiome$mode1$subject), NPLStoolbox::Jakobsen2025$milkMetabolomics$mode1$subject)

# Faeces
newFaeces = NPLStoolbox::Jakobsen2025$faeces

mask = newFaeces$mode1$subject %in% homogenizedSamples
newFaeces$data = newFaeces$data[mask,,]
newFaeces$mode1 = newFaeces$mode1[mask,]
processedFaeces = parafac4microbiome::processDataCube(newFaeces, sparsityThreshold = 0.75, centerMode=1, scaleMode=2)

# Milk
newMilk = NPLStoolbox::Jakobsen2025$milkMicrobiome

mask = newMilk$mode1$subject %in% homogenizedSamples
newMilk$data = newMilk$data[mask,,]
newMilk$mode1 = newMilk$mode1[mask,]
processedMilk = parafac4microbiome::processDataCube(newMilk, sparsityThreshold=0.85, centerMode=1, scaleMode=2)

# Milk metab
newMilkMetab = NPLStoolbox::Jakobsen2025$milkMetabolomics

mask = newMilkMetab$mode1$subject %in% homogenizedSamples
newMilkMetab$data = newMilkMetab$data[mask,,]
newMilkMetab$mode1 = newMilkMetab$mode1[mask,]
processedMilkMetab = parafac4microbiome::processDataCube(newMilkMetab, CLR=FALSE, centerMode=1, scaleMode=2)

datasets = list(processedFaeces$data, processedMilk$data, processedMilkMetab$data)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = setupCMTFdata(datasets, modes, normalize=TRUE)

Y = processedFaeces$mode1$BMI
Ycnt = Y - mean(Y)

pis = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999, 0.9999, 0.99991, 0.99992, 0.99993, 0.99993, 0.99994, 0.99995, 0.99996, 0.99997, 0.99998, 0.99999, 1)
models = list()
for(i in 1:length(pis)){
  pi = pis[i]
  print(pi)
  models[[i]] = CMTFtoolbox::acmtfr_opt(Z, Ycnt, 1, pi=pi, nstart=100, numCores=parallel::detectCores(), method="L-BFGS")
}

saveRDS(models, "CV_BMI_pi_models.RDS")
