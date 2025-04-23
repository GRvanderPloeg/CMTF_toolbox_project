library(tidyverse)
library(parafac4microbiome)
library(NPLStoolbox)
library(CMTFtoolbox)

homogenizedSamples = intersect(intersect(NPLStoolbox::Jakobsen2025$faeces$mode1 %>% filter(!is.na(whz.6m)) %>% select(subject) %>% pull(), NPLStoolbox::Jakobsen2025$milkMicrobiome$mode1 %>% filter(!is.na(whz.6m)) %>% select(subject) %>% pull()), NPLStoolbox::Jakobsen2025$milkMetabolomics$mode1 %>% filter(!is.na(whz.6m)) %>% select(subject) %>% pull())
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

Y = processedFaeces$mode1$whz.6m
Ycnt = Y - mean(Y)
Ynorm = Ycnt / norm(Ycnt, "2")

pi_values = c(0.85, 0.852, 0.854, 0.856, 0.858, 0.86)
for(i in 1:length(pi_values)){
  pi = pi_values[i]
  print(pi)

  result = CMTFtoolbox::ncrossreg(Z, Ynorm, maxNumComponents=10, pi=pi, nstart=10, numCores=parallel::detectCores(), method="L-BFGS", cvFolds=10, abs_tol=1e-6, rel_tol=1e-6, grad_tol=1e-6)
  print(result)
  saveRDS(result, paste0("./CV_small_WHZ_", pi, ".RDS"))
}


