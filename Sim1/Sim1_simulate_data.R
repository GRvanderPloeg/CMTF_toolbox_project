library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)

set.seed(123)

# Settings
delta = 0.05
noiseOnX = 0.30
noiseOnY = 0.05

# Import loadings
inputLoadings = readRDS("./Simulated_data/input_loadings.RDS")

# Sizes of the blocks
numSubjects = nrow(inputLoadings[[1]])
numFeatures1 = nrow(inputLoadings[[2]])
numFeatures2 = nrow(inputLoadings[[4]])
numFeatures3 = nrow(inputLoadings[[6]])
numTimepoints = nrow(inputLoadings[[3]])

# Name matrices
A = inputLoadings[[1]]
B1 = inputLoadings[[2]]
C1 = inputLoadings[[3]]
B2 = inputLoadings[[4]]
C2 = inputLoadings[[5]]
B3 = inputLoadings[[6]]
C3 = inputLoadings[[7]]

# Block X1
X1_term1 = parafac4microbiome::reinflateTensor(A[,1], B1[,1], C1[,1], returnAsTensor=TRUE)
X1_term2 = parafac4microbiome::reinflateTensor(A[,2], B1[,2], C1[,2], returnAsTensor=TRUE)
X1_term3 = parafac4microbiome::reinflateTensor(A[,3], B1[,3], C1[,3], returnAsTensor=TRUE)
X1_term4 = parafac4microbiome::reinflateTensor(A[,4], B1[,4], C1[,4], returnAsTensor=TRUE)
X1_term5 = parafac4microbiome::reinflateTensor(A[,6], B1[,6], C1[,6], returnAsTensor=TRUE)

X1_raw = X1_term1 + X1_term2 + X1_term3 + delta * X1_term4 + X1_term5

# Block X2
X2_term1 = parafac4microbiome::reinflateTensor(A[,1], B2[,1], C2[,1], returnAsTensor=TRUE)
X2_term2 = parafac4microbiome::reinflateTensor(A[,2], B2[,2], C2[,2], returnAsTensor=TRUE)
X2_term3 = parafac4microbiome::reinflateTensor(A[,3], B2[,3], C2[,3], returnAsTensor=TRUE)
X2_term4 = parafac4microbiome::reinflateTensor(A[,5], B2[,5], C2[,5], returnAsTensor=TRUE)
X2_term5 = parafac4microbiome::reinflateTensor(A[,7], B2[,7], C2[,7], returnAsTensor=TRUE)

X2_raw = X2_term1 + X2_term2 + X2_term3 + X2_term4 + X2_term5

# Block X3
X3_term1 = parafac4microbiome::reinflateTensor(A[,1], B3[,1], C3[,1], returnAsTensor=TRUE)
X3_term2 = parafac4microbiome::reinflateTensor(A[,2], B3[,2], C3[,2], returnAsTensor=TRUE)
X3_term3 = parafac4microbiome::reinflateTensor(A[,4], B3[,4], C3[,4], returnAsTensor=TRUE)
X3_term4 = parafac4microbiome::reinflateTensor(A[,5], B3[,5], C3[,5], returnAsTensor=TRUE)
X3_term5 = parafac4microbiome::reinflateTensor(A[,8], B3[,8], C3[,8], returnAsTensor=TRUE)

X3_raw = X3_term1 + X3_term2 + delta * X3_term3 + X3_term4 + X3_term5

# Create noise
noise1 = as.tensor(array(rnorm(numSubjects*numFeatures1*numTimepoints), c(numSubjects, numFeatures1, numTimepoints)))
noise2 = as.tensor(array(rnorm(numSubjects*numFeatures2*numTimepoints), c(numSubjects, numFeatures2, numTimepoints)))
noise3 = as.tensor(array(rnorm(numSubjects*numFeatures3*numTimepoints), c(numSubjects, numFeatures3, numTimepoints)))

# Make the noise the right frobenius norm
noise1_scaled = fnorm(X1_raw) / (1/noiseOnX) * (noise1 / fnorm(noise1))
noise2_scaled = fnorm(X2_raw) / (1/noiseOnX) * (noise2 / fnorm(noise2))
noise3_scaled = fnorm(X3_raw) / (1/noiseOnX) * (noise3 / fnorm(noise3))

# Add noise to X blocks
X1_final = X1_raw + noise1_scaled
X2_final = X2_raw + noise2_scaled
X3_final = X3_raw + noise3_scaled

# Generate Y
Y = A[,4]
noiseY = as.matrix(rnorm(numSubjects))
noiseY = noiseY - mean(noiseY)
noiseY_scaled = (norm(Y, "2") / (1/noiseOnY)) * (noiseY / norm(noiseY, "2"))
Y_with_noise = Y + noiseY_scaled
Y_final = Y_with_noise / norm(Y_with_noise, "2")
Y_final = as.tensor(Y_final)

# Save sim
saveRDS(X1_final@data, "./Sim1/Sim1_X1.RDS")
saveRDS(X2_final@data, "./Sim1/Sim1_X2.RDS")
saveRDS(X3_final@data, "./Sim1/Sim1_X3.RDS")
saveRDS(Y_final@data, "./Sim1/Sim1_Y.RDS")
