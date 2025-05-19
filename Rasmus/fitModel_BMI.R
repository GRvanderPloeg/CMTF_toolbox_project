library(tidyverse)
library(parafac4microbiome)
library(NPLStoolbox)
library(CMTFtoolbox)

homogenizedSamples = intersect(intersect(NPLStoolbox::Jakobsen2025$faeces$mode1 %>% filter(!is.na(BMI)) %>% select(subject) %>% pull(), NPLStoolbox::Jakobsen2025$milkMicrobiome$mode1 %>% filter(!is.na(BMI)) %>% select(subject) %>% pull()), NPLStoolbox::Jakobsen2025$milkMetabolomics$mode1 %>% filter(!is.na(BMI)) %>% select(subject) %>% pull())
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
Ynorm = Ycnt / norm(Ycnt, "2")

model = CMTFtoolbox::acmtfr_opt(Z, Ynorm, numComponents=1, pi=0.90, nstart=100, numCores=parallel::detectCores())
saveRDS(model, "./ACMTFR_model_bmi90.RDS")

model = CMTFtoolbox::acmtfr_opt(Z, Ynorm, numComponents=1, pi=0.85, nstart=100, numCores=parallel::detectCores())
saveRDS(model, "./ACMTFR_model_bmi85.RDS")

model = CMTFtoolbox::acmtfr_opt(Z, Ynorm, numComponents=1, pi=0.80, nstart=100, numCores=parallel::detectCores())
saveRDS(model, "./ACMTFR_model_bmi80.RDS")

model = CMTFtoolbox::acmtfr_opt(Z, Ynorm, numComponents=1, pi=0.75, nstart=100, numCores=parallel::detectCores())
saveRDS(model, "./ACMTFR_model_bmi75.RDS")
