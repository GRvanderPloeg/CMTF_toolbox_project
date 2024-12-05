library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)

library(parallel)
library(doParallel)
library(foreach)
set.seed(123)

# Sizes of the blocks
numSubjects = 108
numFeatures = 100
numTimepoints = 10

# Settings
Ysize = 1
noiseOnX = 0.25 # Noise percentage on each X block
noiseOnY = 0.05 # Noise percentage on Y
w_global = 99
w_local = 1
w_distinct1 = 5
w_distinct2 = 5
w_distinct3 = 5

# Simulation
r = array(rnorm(numSubjects*5), c(numSubjects, 5))
r = sweep(r, 2, colMeans(r), FUN="-")
U = svd(r)$u # creates orthogonal score vectors with norm 1 mean 0

a_global = U[,1]
a_local = U[,2]
a_distinct1 = U[,3]
a_distinct2 = U[,4]
a_distinct3 = U[,5]

b_global1 = array(0L, c(numFeatures, 1))
b_global1[1:4] = 1
b_global1 = b_global1 - mean(b_global1)
b_global1 = b_global1 / norm(b_global1, "2")

b_global2 = array(0L, c(numFeatures, 1))
b_global2[5:8] = 1
b_global2 = b_global2 - mean(b_global2)
b_global2 = b_global2 / norm(b_global2, "2")

b_global3 = array(0L, c(numFeatures, 1))
b_global3[9:12] = 1
b_global3 = b_global3 - mean(b_global3)
b_global3 = b_global3 / norm(b_global3, "2")

b_local1 = array(0L, c(numFeatures, 1))
b_local1[21:26] = 1
b_local1 = b_local1 - mean(b_local1)
b_local1 = b_local1 / norm(b_local1, "2")

b_local2 = array(0L, c(numFeatures, 1))
b_local2[27:32] = 1
b_local2 = b_local2 - mean(b_local2)
b_local2 = b_local2 / norm(b_local2, "2")

# Block 3 does not load onto the local component, as simplification

b_distinct1 = array(0L, c(numFeatures, 1))
b_distinct1[41:52] = 1
b_distinct1 = b_distinct1 - mean(b_distinct1)
b_distinct1 = b_distinct1 / norm(b_distinct1, "2")

b_distinct2 = array(0L, c(numFeatures, 1))
b_distinct2[60:72] = 1
b_distinct2 = b_distinct2 - mean(b_distinct2)
b_distinct2 = b_distinct2 / norm(b_distinct2, "2")

b_distinct3 = array(0L, c(numFeatures, 1))
b_distinct3[80:92] = 1
b_distinct3 = b_distinct3 - mean(b_distinct3)
b_distinct3 = b_distinct3 / norm(b_distinct3, "2")

c_global1 = t(exp(1:numTimepoints))
c_global1 = c_global1 - mean(c_global1)
c_global1 = c_global1 / norm(c_global1, "2")
c_global1 = t(c_global1)

c_global2 = t(sin(1:numTimepoints))
c_global2 = c_global2 - mean(c_global2)
c_global2 = c_global2 / norm(c_global2, "2")
c_global2 = t(c_global2)

c_global3 = t(cos(1:numTimepoints))
c_global3 = c_global3 - mean(c_global3)
c_global3 = c_global3 / norm(c_global3, "2")
c_global3 = t(c_global3)

c_local1 = t(sqrt(1:numTimepoints))
c_local1 = c_local1 - mean(c_local1)
c_local1 = c_local1 / norm(c_local1, "2")
c_local1 = t(c_local1)

c_local2 = t((1:numTimepoints)^2)
c_local2 = c_local2 - mean(c_local2)
c_local2 = c_local2 / norm(c_local2, "2")
c_local2 = t(c_local2)

c_distinct1 = t(1:numTimepoints)
c_distinct1 = c_distinct1 - mean(c_distinct1)
c_distinct1 = c_distinct1 / norm(c_distinct1, "2")
c_distinct1 = t(c_distinct1)

c_distinct2 = t(1/(1:numTimepoints))
c_distinct2 = c_distinct2 - mean(c_distinct2)
c_distinct2 = c_distinct2 / norm(c_distinct2, "2")
c_distinct2 = t(c_distinct2)

c_distinct3 = t(log(1:numTimepoints))
c_distinct3 = c_distinct3 - mean(c_distinct3)
c_distinct3 = c_distinct3 / norm(c_distinct3, "2")
c_distinct3 = t(c_distinct3)

# Block X1
X1_term1 = parafac4microbiome::reinflateTensor(a_global, b_global1, c_global1, returnAsTensor=TRUE) # rTensor::fnorm is 1
X1_term2 = parafac4microbiome::reinflateTensor(a_local, b_local1, c_local1, returnAsTensor=TRUE) # rTensor::fnorm is 1
X1_term3 = parafac4microbiome::reinflateTensor(a_distinct1, b_distinct1, c_distinct1, returnAsTensor=TRUE) # rTensor::fnorm is 1

X1_raw = w_global * X1_term1 + w_local * X1_term2 + w_distinct1 * X1_term3

# Block X2
X2_term1 = parafac4microbiome::reinflateTensor(a_global, b_global2, c_global2, returnAsTensor=TRUE) # rTensor::fnorm is 1
X2_term2 = parafac4microbiome::reinflateTensor(a_local, b_local2, c_local2, returnAsTensor=TRUE) # rTensor::fnorm is 1
X2_term3 = parafac4microbiome::reinflateTensor(a_distinct2, b_distinct2, c_distinct2, returnAsTensor=TRUE) # rTensor::fnorm is 1

X2_raw = w_global * X2_term1 + w_local * X2_term2 + w_distinct2 * X2_term3

# Block X3
X3_term1 = parafac4microbiome::reinflateTensor(a_global, b_global3, c_global3, returnAsTensor=TRUE) # rTensor::fnorm is 1
# there is no local component
X3_term3 = parafac4microbiome::reinflateTensor(a_distinct3, b_distinct3, c_distinct3, returnAsTensor=TRUE) # rTensor::fnorm is 1

X3_raw = w_global * X3_term1 + w_distinct3 * X3_term3

# Create noise
noise1 = as.tensor(array(rnorm(numSubjects*numFeatures*numTimepoints), c(numSubjects, numFeatures, numTimepoints)))
noise2 = as.tensor(array(rnorm(numSubjects*numFeatures*numTimepoints), c(numSubjects, numFeatures, numTimepoints)))
noise3 = as.tensor(array(rnorm(numSubjects*numFeatures*numTimepoints), c(numSubjects, numFeatures, numTimepoints)))

# Center the noise
noise1 = parafac4microbiome::multiwayCenter(noise1)
noise2 = parafac4microbiome::multiwayCenter(noise2)
noise3 = parafac4microbiome::multiwayCenter(noise3)

# Scale the noise
noise1 = parafac4microbiome::multiwayScale(noise1)
noise2 = parafac4microbiome::multiwayScale(noise2)
noise3 = parafac4microbiome::multiwayScale(noise3)

# Put back into Tensor
noise1 = as.tensor(noise1)
noise2 = as.tensor(noise2)
noise3 = as.tensor(noise3)

# Make the noise the right frobenius norm
noise1 = fnorm(X1_raw) / (1/noiseOnX) * (noise1 / fnorm(noise1))
noise2 = fnorm(X2_raw) / (1/noiseOnX) * (noise2 / fnorm(noise2))
noise3 = fnorm(X3_raw) / (1/noiseOnX) * (noise3 / fnorm(noise3))

# Add noise to X blocks
X1_noise = X1_raw + noise1
X2_noise = X2_raw + noise2
X3_noise = X3_raw + noise3

# Center
X1_cnt = parafac4microbiome::multiwayCenter(X1_noise)
X2_cnt = parafac4microbiome::multiwayCenter(X2_noise)
X3_cnt = parafac4microbiome::multiwayCenter(X3_noise)

# Scale
X1_cnt_scl = parafac4microbiome::multiwayScale(X1_cnt)
X2_cnt_scl = parafac4microbiome::multiwayScale(X2_cnt)
X3_cnt_scl = parafac4microbiome::multiwayScale(X3_cnt)

# Define final version of X
X1_final = X1_cnt_scl
X2_final = X2_cnt_scl
X3_final = X3_cnt_scl

# Prepare Y
Y = a_local
noiseY = rnorm(numSubjects)
noiseY = noiseY - mean(noiseY)
noiseY = (norm(Y, "2") / (1/noiseOnY)) * (noiseY / norm(noiseY, "2"))
Ynoise = Y + noiseY
Ynorm = Ynoise / norm(Ynoise, "2")
Y_final = Ysize * Ynorm

# Convert to Tensor
Y_final = as.tensor(Y_final)

# Prepare Z
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=TRUE)

# ACMTFR
betas = c(5e-5, 1e-5, 5e-4, 1e-4, 5e-3, 1e-3, 5e-2, 1e-2, 5e-1, 1e-1)
pis   = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

simSettings = cbind(rep(betas, each=length(pis)), rep(pis, length(betas)))

cl = parallel::makeCluster(25)
doParallel::registerDoParallel(cl)
models = foreach::foreach(i=1:nrow(simSettings)) %dopar% {
  betaValue = simSettings[i,1]
  piValue = simSettings[i,2]
  model=CMTFtoolbox::acmtfr_opt(Z,matrix(Y_final@data),numComponents=5,initialization="nvec",beta=rep(betaValue,3),pi=piValue,abs_tol=1e-8,rel_tol=1e-8,nstart=1)
}
parallel::stopCluster(cl)

saveRDS(models, "./20241121_ACMTFR_simulation.RDS")


