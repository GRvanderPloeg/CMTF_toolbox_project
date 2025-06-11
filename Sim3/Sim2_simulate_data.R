library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)

set.seed(123)

# Sizes of the blocks
numSubjects = 100
numFeatures1 = 120
numFeatures2 = 150
numFeatures3 = 90
numTimepoints = 5

# Settings
delta = 0
# gamma = 0.001 # varies

noiseOnX = 0.50
noiseOnY = 0

# Simulation

# initialize subject loadings
r = array(rnorm(numSubjects*7), c(numSubjects, 7))
U = CMTFtoolbox::removeTwoNormCol(r)

a_global1 = U[,1]
a_global2 = U[,2]
a_local1 = U[,3]
a_local2 = U[,4]
a_distinct1 = U[,5]
a_distinct2 = U[,6]
a_distinct3 = U[,7]
scores = U

# initialize feature loadings 1
r = array(rnorm(numFeatures1*7), c(numFeatures1, 7))
U = CMTFtoolbox::removeTwoNormCol(r)

b1_global1 = U[,1]
b1_global2 = U[,2]
b1_local1 = U[,3]
b1_local2 = U[,4]
b1_distinct1 = U[,5]
b1_distinct2 = U[,6]
b1_distinct3 = U[,7]

loadings1 = U

# initialize feature loadings 2
r = array(rnorm(numFeatures2*7), c(numFeatures2, 7))
U = CMTFtoolbox::removeTwoNormCol(r)

b2_global1 = U[,1]
b2_global2 = U[,2]
b2_local1 = U[,3]
b2_local2 = U[,4]
b2_distinct1 = U[,5]
b2_distinct2 = U[,6]
b2_distinct3 = U[,7]

loadings2 = U

# initialize feature loadings 3
r = array(rnorm(numFeatures3*7), c(numFeatures3, 7))
U = CMTFtoolbox::removeTwoNormCol(r)

b3_global1 = U[,1]
b3_global2 = U[,2]
b3_local1 = U[,3]
b3_local2 = U[,4]
b3_distinct1 = U[,5]
b3_distinct2 = U[,6]
b3_distinct3 = U[,7]

loadings3 = U

# initialize time loadings 1
r = array(rnorm(numTimepoints*7), c(numTimepoints, 7))
U = CMTFtoolbox::removeTwoNormCol(r)

c1_global1 = U[,1]
c1_global2 = U[,2]
c1_local1 = U[,3]
c1_local2 = U[,4]
c1_distinct1 = U[,5]
c1_distinct2 = U[,6]
c1_distinct3 = U[,7]

timeLoadings1 = U

# initialize time loadings 2
r = array(rnorm(numTimepoints*7), c(numTimepoints, 7))
U = CMTFtoolbox::removeTwoNormCol(r)

c2_global1 = U[,1]
c2_global2 = U[,2]
c2_local1 = U[,3]
c2_local2 = U[,4]
c2_distinct1 = U[,5]
c2_distinct2 = U[,6]
c2_distinct3 = U[,7]

timeLoadings2 = U

# initialize time loadings 3
r = array(rnorm(numTimepoints*7), c(numTimepoints, 7))
U = CMTFtoolbox::removeTwoNormCol(r)

c3_global1 = U[,1]
c3_global2 = U[,2]
c3_local1 = U[,3]
c3_local2 = U[,4]
c3_distinct1 = U[,5]
c3_distinct2 = U[,6]
c3_distinct3 = U[,7]

timeLoadings3 = U

# Block X1
X1_term1 = parafac4microbiome::reinflateTensor(a_global1, b1_global1, c1_global1, returnAsTensor=TRUE)
X1_term2 = parafac4microbiome::reinflateTensor(a_global2, b1_global2, c1_global2, returnAsTensor=TRUE)
X1_term3 = parafac4microbiome::reinflateTensor(a_local1, b1_local1, c1_local1, returnAsTensor=TRUE)
X1_term4 = parafac4microbiome::reinflateTensor(a_local2, b1_local2, c1_local2, returnAsTensor=TRUE)
X1_term5 = parafac4microbiome::reinflateTensor(a_distinct1, b1_distinct1, c1_distinct1, returnAsTensor=TRUE)
X1_term6 = parafac4microbiome::reinflateTensor(a_distinct2, b1_distinct2, c1_distinct2, returnAsTensor=TRUE)
X1_term7 = parafac4microbiome::reinflateTensor(a_distinct3, b1_distinct3, c1_distinct3, returnAsTensor=TRUE)

# Block X2
X2_term1 = parafac4microbiome::reinflateTensor(a_global1, b2_global1, c2_global1, returnAsTensor=TRUE)
X2_term2 = parafac4microbiome::reinflateTensor(a_global2, b2_global2, c2_global2, returnAsTensor=TRUE)
X2_term3 = parafac4microbiome::reinflateTensor(a_local1, b2_local1, c2_local1, returnAsTensor=TRUE)
X2_term4 = parafac4microbiome::reinflateTensor(a_local2, b2_local2, c2_local2, returnAsTensor=TRUE)
X2_term5 = parafac4microbiome::reinflateTensor(a_distinct1, b2_distinct1, c2_distinct1, returnAsTensor=TRUE)
X2_term6 = parafac4microbiome::reinflateTensor(a_distinct2, b2_distinct2, c2_distinct2, returnAsTensor=TRUE)
X2_term7 = parafac4microbiome::reinflateTensor(a_distinct3, b2_distinct3, c2_distinct3, returnAsTensor=TRUE)

# Block X3
X3_term1 = parafac4microbiome::reinflateTensor(a_global1, b3_global1, c3_global1, returnAsTensor=TRUE)
X3_term2 = parafac4microbiome::reinflateTensor(a_global2, b3_global2, c3_global2, returnAsTensor=TRUE)
X3_term3 = parafac4microbiome::reinflateTensor(a_local1, b3_local1, c3_local1, returnAsTensor=TRUE)
X3_term4 = parafac4microbiome::reinflateTensor(a_local2, b3_local2, c3_local2, returnAsTensor=TRUE)
X3_term5 = parafac4microbiome::reinflateTensor(a_distinct1, b3_distinct1, c3_distinct1, returnAsTensor=TRUE)
X3_term6 = parafac4microbiome::reinflateTensor(a_distinct2, b3_distinct2, c3_distinct2, returnAsTensor=TRUE)
X3_term7 = parafac4microbiome::reinflateTensor(a_distinct3, b3_distinct3, c3_distinct3, returnAsTensor=TRUE)

# Create noise
noise1 = as.tensor(array(rnorm(numSubjects*numFeatures1*numTimepoints), c(numSubjects, numFeatures1, numTimepoints)))
noise2 = as.tensor(array(rnorm(numSubjects*numFeatures2*numTimepoints), c(numSubjects, numFeatures2, numTimepoints)))
noise3 = as.tensor(array(rnorm(numSubjects*numFeatures3*numTimepoints), c(numSubjects, numFeatures3, numTimepoints)))

# Create the blocks
gammas = c(0.01, 0.1, 0.25, 0.50, 0.75, 1)

for(i in 1:length(gammas)){
  gamma = gammas[i]

  X1_raw = X1_term1 + X1_term2 + X1_term3 + gamma * X1_term4 + X1_term5 + delta * X1_term6 + delta * X1_term7
  X2_raw = X2_term1 + X2_term2 + X2_term3 + gamma * X2_term4 + delta * X2_term5 + X2_term6 + delta * X2_term7
  X3_raw = X3_term1 + X3_term2 + delta * X3_term3 + delta * X3_term4 + delta * X3_term5 + delta * X3_term6 + X3_term7

  # Make the noise the right frobenius norm
  noise1_scaled = fnorm(X1_raw) / (1/noiseOnX) * (noise1 / fnorm(noise1))
  noise2_scaled = fnorm(X2_raw) / (1/noiseOnX) * (noise2 / fnorm(noise2))
  noise3_scaled = fnorm(X3_raw) / (1/noiseOnX) * (noise3 / fnorm(noise3))

  # Add noise to X blocks
  X1_final = X1_raw + noise1_scaled
  X2_final = X2_raw + noise2_scaled
  X3_final = X3_raw + noise3_scaled

  saveRDS(X1_final@data, paste0("./Sim2/Sim2_X1_", gamma, ".RDS"))
  saveRDS(X2_final@data, paste0("./Sim2/Sim2_X2_", gamma, ".RDS"))
  saveRDS(X3_final@data, paste0("./Sim2/Sim2_X3_", gamma, ".RDS"))
}

# Generate Y
Y_final = as.tensor(a_local2)
# noiseY = rnorm(numSubjects)
# noiseY = noiseY - mean(noiseY)
# noiseY = (norm(Y, "2") / (1/noiseOnY)) * (noiseY / norm(noiseY, "2"))
# Ynoise = Y + noiseY
# Y_final = Ynoise / norm(Ynoise, "2")
# Y_final = as.tensor(Y_final)

# Save Y
saveRDS(Y_final@data, "./Sim2/Sim2_Y.RDS")

# Save loadings
allLoadings = list(scores, loadings1, loadings2, loadings3, timeLoadings1, timeLoadings2, timeLoadings3)
saveRDS(allLoadings, "./Sim2/Sim2_input_loadings.RDS")
