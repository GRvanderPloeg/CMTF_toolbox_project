library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)

set.seed(123)

# Sizes of the blocks
numSubjects = 50
numFeatures = 30
numTimepoints = 40

# Settings
Ysize = 1
noiseOnX = 0.25 # Noise percentage on each X block
#noiseOnY = 0.01 # Set lower down
w_global1 = 1
w_global2 = 1
w_local1 = 1
w_local2 = 1
w_distinct1 = 1
w_distinct2 = 1
w_distinct3 = 1
w_other = 0.05

rho1 = 0
rho2 = 0
rho3 = 1
rho4 = 0

# Simulation

# initialize subject loadings
r = array(rnorm(numSubjects*7), c(numSubjects, 7))
U = CMTFtoolbox::removeTwoNormCol(r)
# r = sweep(r, 2, colMeans(r), FUN="-")
# U = svd(r)$u

a_global1 = U[,1]
a_global2 = U[,2]
a_local1 = U[,3]
a_local2 = U[,4]
a_distinct1 = U[,5]
a_distinct2 = U[,6]
a_distinct3 = U[,7]
scores = U

# initialize feature loadings
r = array(rnorm(numFeatures*21), c(numFeatures, 21))
U = CMTFtoolbox::removeTwoNormCol(r)
# r = sweep(r, 2, colMeans(r), FUN="-")
# U = svd(r)$u

b_global1_X1 = U[,1]
b_global1_X2 = U[,2]
b_global1_X3 = U[,3]

b_global2_X1 = U[,4]
b_global2_X2 = U[,5]
b_global2_X3 = U[,6]

b_local1_X1 = U[,7]
b_local1_X2 = U[,8]
b_local1_X3 = U[,9]

b_local2_X1 = U[,10]
b_local2_X2 = U[,11]
b_local2_X3 = U[,12]

b_distinct1_X1 = U[,13]
b_distinct1_X2 = U[,14]
b_distinct1_X3 = U[,15]

b_distinct2_X1 = U[,16]
b_distinct2_X2 = U[,17]
b_distinct2_X3 = U[,18]

b_distinct3_X1 = U[,19]
b_distinct3_X2 = U[,20]
b_distinct3_X3 = U[,21]

loadings = U

# initialize time loadings
r = array(rnorm(numTimepoints*21), c(numTimepoints, 21))
U = CMTFtoolbox::removeTwoNormCol(r)
# r = sweep(r, 2, colMeans(r), FUN="-")
# U = svd(r)$u

c_global1_X1 = U[,1]
c_global1_X2 = U[,2]
c_global1_X3 = U[,3]

c_global2_X1 = U[,4]
c_global2_X2 = U[,5]
c_global2_X3 = U[,6]

c_local1_X1 = U[,7]
c_local1_X2 = U[,8]
c_local1_X3 = U[,9]

c_local2_X1 = U[,10]
c_local2_X2 = U[,11]
c_local2_X3 = U[,12]

c_distinct1_X1 = U[,13]
c_distinct1_X2 = U[,14]
c_distinct1_X3 = U[,15]

c_distinct2_X1 = U[,16]
c_distinct2_X2 = U[,17]
c_distinct2_X3 = U[,18]

c_distinct3_X1 = U[,19]
c_distinct3_X2 = U[,20]
c_distinct3_X3 = U[,21]

timeLoadings = U

# Block X1
X1_term1 = parafac4microbiome::reinflateTensor(a_global1, b_global1_X1, c_global1_X1, returnAsTensor=TRUE) # rTensor::fnorm is 1
X1_term2 = parafac4microbiome::reinflateTensor(a_global2, b_global2_X1, c_global2_X1, returnAsTensor=TRUE) # rTensor::fnorm is 1
X1_term3 = parafac4microbiome::reinflateTensor(a_local1, b_local1_X1, c_local1_X1, returnAsTensor=TRUE) # rTensor::fnorm is 1
X1_term4 = parafac4microbiome::reinflateTensor(a_local2, b_local2_X1, c_local2_X1, returnAsTensor=TRUE) # rTensor::fnorm is 1
X1_term5 = parafac4microbiome::reinflateTensor(a_distinct1, b_distinct1_X1, c_distinct1_X1, returnAsTensor=TRUE) # rTensor::fnorm is 1
X1_term6 = parafac4microbiome::reinflateTensor(a_distinct2, b_distinct2_X1, c_distinct2_X1, returnAsTensor=TRUE) # rTensor::fnorm is 1
X1_term7 = parafac4microbiome::reinflateTensor(a_distinct3, b_distinct3_X1, c_distinct3_X1, returnAsTensor=TRUE) # rTensor::fnorm is 1

X1_raw = w_global1 * X1_term1 + w_global2 * X1_term2 + w_local1 * X1_term3 + w_local2 * X1_term4 + w_distinct1 * X1_term5 + w_other * X1_term6 + w_other * X1_term7

# Block X2
X2_term1 = parafac4microbiome::reinflateTensor(a_global1, b_global1_X2, c_global1_X2, returnAsTensor=TRUE) # rTensor::fnorm is 1
X2_term2 = parafac4microbiome::reinflateTensor(a_global2, b_global2_X2, c_global2_X2, returnAsTensor=TRUE) # rTensor::fnorm is 1
X2_term3 = parafac4microbiome::reinflateTensor(a_local1, b_local1_X2, c_local1_X2, returnAsTensor=TRUE) # rTensor::fnorm is 1
X2_term4 = parafac4microbiome::reinflateTensor(a_local2, b_local2_X2, c_local2_X2, returnAsTensor=TRUE) # rTensor::fnorm is 1
X2_term5 = parafac4microbiome::reinflateTensor(a_distinct1, b_distinct1_X2, c_distinct1_X2, returnAsTensor=TRUE) # rTensor::fnorm is 1
X2_term6 = parafac4microbiome::reinflateTensor(a_distinct2, b_distinct2_X2, c_distinct2_X2, returnAsTensor=TRUE) # rTensor::fnorm is 1
X2_term7 = parafac4microbiome::reinflateTensor(a_distinct3, b_distinct3_X2, c_distinct3_X2, returnAsTensor=TRUE) # rTensor::fnorm is 1

X2_raw = w_global1 * X2_term1 + w_global2 * X2_term2 + w_local1 * X2_term3 + w_local2 * X2_term4 + w_other * X2_term5 +  w_distinct2 * X2_term6 + w_other * X2_term7

# Block X3
X3_term1 = parafac4microbiome::reinflateTensor(a_global1, b_global1_X3, c_global1_X3, returnAsTensor=TRUE) # rTensor::fnorm is 1
X3_term2 = parafac4microbiome::reinflateTensor(a_global2, b_global2_X3, c_global2_X3, returnAsTensor=TRUE) # rTensor::fnorm is 1
X3_term3 = parafac4microbiome::reinflateTensor(a_local1, b_local1_X3, c_local1_X3, returnAsTensor=TRUE) # rTensor::fnorm is 1
X3_term4 = parafac4microbiome::reinflateTensor(a_local2, b_local2_X3, c_local2_X3, returnAsTensor=TRUE) # rTensor::fnorm is 1
X3_term5 = parafac4microbiome::reinflateTensor(a_distinct1, b_distinct1_X3, c_distinct1_X3, returnAsTensor=TRUE) # rTensor::fnorm is 1
X3_term6 = parafac4microbiome::reinflateTensor(a_distinct2, b_distinct2_X3, c_distinct2_X3, returnAsTensor=TRUE) # rTensor::fnorm is 1
X3_term7 = parafac4microbiome::reinflateTensor(a_distinct3, b_distinct3_X3, c_distinct3_X3, returnAsTensor=TRUE) # rTensor::fnorm is 1

X3_raw = w_global1 * X3_term1 + w_global2 * X3_term2 + w_other * X3_term3 + w_other * X3_term4 + w_other * X3_term5 + w_other * X3_term6 + w_distinct3 * X3_term7

# Create noise
noise1 = as.tensor(array(rnorm(numSubjects*numFeatures*numTimepoints), c(numSubjects, numFeatures, numTimepoints)))
noise2 = as.tensor(array(rnorm(numSubjects*numFeatures*numTimepoints), c(numSubjects, numFeatures, numTimepoints)))
noise3 = as.tensor(array(rnorm(numSubjects*numFeatures*numTimepoints), c(numSubjects, numFeatures, numTimepoints)))

# Center the noise
# noise1 = parafac4microbiome::multiwayCenter(noise1)
# noise2 = parafac4microbiome::multiwayCenter(noise2)
# noise3 = parafac4microbiome::multiwayCenter(noise3)
#
# # Scale the noise
# noise1 = parafac4microbiome::multiwayScale(noise1)
# noise2 = parafac4microbiome::multiwayScale(noise2)
# noise3 = parafac4microbiome::multiwayScale(noise3)

# # Put back into Tensor
# noise1 = as.tensor(noise1)
# noise2 = as.tensor(noise2)
# noise3 = as.tensor(noise3)

# Make the noise the right frobenius norm
noise1 = fnorm(X1_raw) / (1/noiseOnX) * (noise1 / fnorm(noise1))
noise2 = fnorm(X2_raw) / (1/noiseOnX) * (noise2 / fnorm(noise2))
noise3 = fnorm(X3_raw) / (1/noiseOnX) * (noise3 / fnorm(noise3))

# Add noise to X blocks
X1_noise = X1_raw + noise1
X2_noise = X2_raw + noise2
X3_noise = X3_raw + noise3

# # Center
# X1_cnt = parafac4microbiome::multiwayCenter(X1_noise)
# X2_cnt = parafac4microbiome::multiwayCenter(X2_noise)
# X3_cnt = parafac4microbiome::multiwayCenter(X3_noise)
#
# # Scale
# X1_cnt_scl = parafac4microbiome::multiwayScale(X1_cnt)
# X2_cnt_scl = parafac4microbiome::multiwayScale(X2_cnt)
# X3_cnt_scl = parafac4microbiome::multiwayScale(X3_cnt)
#
# # Define final version of X
# X1_final = X1_cnt_scl
# X2_final = X2_cnt_scl
# X3_final = X3_cnt_scl

X1_final = X1_noise@data
X2_final = X2_noise@data
X3_final = X3_noise@data

# Generate Y
Y = rho1 * a_global1 + rho2 * a_global2 + rho3 * a_local1 + rho4 * a_local2
Ynoise = Y
Ynorm = Ynoise / norm(Ynoise, "2")
Y_final = Ysize * Ynorm
Y_final = as.tensor(Y_final)

saveRDS(Y_final, paste0("./F_test/Sim2b_Y_Y_inside.RDS"))

# Save sim
saveRDS(X1_final, "./F_test/Sim2b_X1_Y_inside.RDS")
saveRDS(X2_final, "./F_test/Sim2b_X2_Y_inside.RDS")
saveRDS(X3_final, "./F_test/Sim2b_X3_Y_inside.RDS")

allLoadings = list(scores, loadings, timeLoadings)
saveRDS(allLoadings, "./F_test/Sim2b_input_loadings_Y_inside.RDS")

parameters = c(w_global1, w_global2, w_local1, w_local2, w_distinct1, w_distinct2, w_distinct3, rho1, rho2, rho3, rho4, noiseOnX, Ysize)
saveRDS(parameters, "./F_test/Sim2b_parameters_Y_inside.RDS")
