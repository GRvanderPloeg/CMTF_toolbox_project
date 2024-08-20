library(CMTFtoolbox)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rTensor)
library(parallel)

set.seed(123)

FMS_to_simulation = function(inputFac, Fac, modes){
  numComponents = ncol(Fac[[1]])
  numModes = max(unlist(modes))
  numDatasets = length(modes)
  normalizedInputFac = normalizeFac(inputFac, modes)
  normalizedFac = normalizeFac(Fac, modes)

  # Find which components match between the model and the input
  matchingComponents = apply(abs(cor(normalizedInputFac$Fac[[1]], normalizedFac$Fac[[1]])), 2, which.max)

  # Calculate FMS
  FMS_result = 1:numComponents
  for(i in 1:numComponents){
    FMS = 1
    ksi = 0
    ksi_hat = 0
    matchingComponent = matchingComponents[i]

    for(j in 1:numModes){
      vect1 = normalizedInputFac$Fac[[j]][,matchingComponent]
      vect2 = normalizedFac$Fac[[j]][,i]
      FMS = FMS * (t(vect1) %*% vect2)
    }

    for(j in 1:numDatasets){
      ksi = ksi + normalizedInputFac$normsPerDataset[j,matchingComponent]
      ksi_hat = ksi_hat + normalizedFac$normsPerDataset[j,i]
    }
    ksi_term = ((abs(ksi - ksi_hat)) / max(ksi, ksi_hat))
    FMS_result[i] = (1 - ksi_term) * abs(FMS)
  }
  return(min(FMS_result))
}

R = 3
I = 100
J = 100
K = 100
L = 100
M = 100
modes = list(c(1,2,3), c(1,4,5))

# Generate loading matrices, normalize each vector to 1
A = array(rnorm(I*R), c(I, R))  # shared subject mode
Anorm = sweep(A, 2, apply(A, 2, function(x){norm(as.matrix(x), "F")}), FUN="/")
B = array(rnorm(J*R), c(J, R))  # distinct feature mode of X1
Bnorm = sweep(B, 2, apply(B, 2, function(x){norm(as.matrix(x), "F")}), FUN="/")
C = array(rnorm(K*R), c(K, R))  # distinct condition mode of X1
Cnorm = sweep(C, 2, apply(C, 2, function(x){norm(as.matrix(x), "F")}), FUN="/")
D = array(rnorm(L*R), c(L, R))  # distinct feature mode of X2
Dnorm = sweep(D, 2, apply(D, 2, function(x){norm(as.matrix(x), "F")}), FUN="/")
E = array(rnorm(M*R), c(M, R))  # distinct condition mode of X2
Enorm = sweep(E, 2, apply(E, 2, function(x){norm(as.matrix(x), "F")}), FUN="/")

# Generate datasets without noise
df1 = as.tensor(reinflateTensor(A, B, C))
df2 = as.tensor(reinflateTensor(A, D, E))

df1norm = as.tensor(reinflateTensor(Anorm, Bnorm, Cnorm))
df2norm = as.tensor(reinflateTensor(Anorm, Dnorm, Enorm))

# Generate noise
df1_noise = as.tensor(array(rnorm(I*J*K), c(I,J,K)))
df2_noise = as.tensor(array(rnorm(I*L*M), c(I,L,M)))

# Create final datasets by adding the noise of different sizes to the generated data
X1_010 = df1 + 0.10 * df1_noise * fnorm(df1) / fnorm(df1_noise)
X1_025 = df1 + 0.25 * df1_noise * fnorm(df1) / fnorm(df1_noise)
X1_035 = df1 + 0.35 * df1_noise * fnorm(df1) / fnorm(df1_noise)

X2_010 = df2 + 0.10 * df2_noise * fnorm(df2) / fnorm(df2_noise)
X2_025 = df2 + 0.25 * df2_noise * fnorm(df2) / fnorm(df2_noise)
X2_035 = df2 + 0.35 * df2_noise * fnorm(df2) / fnorm(df2_noise)

X1_010_norm = df1norm + 0.10 * df1_noise * fnorm(df1norm) / fnorm(df1_noise)
X1_025_norm = df1norm + 0.25 * df1_noise * fnorm(df1norm) / fnorm(df1_noise)
X1_035_norm = df1norm + 0.35 * df1_noise * fnorm(df1norm) / fnorm(df1_noise)

X2_010_norm = df2norm + 0.10 * df2_noise * fnorm(df2norm) / fnorm(df2_noise)
X2_025_norm = df2norm + 0.25 * df2_noise * fnorm(df2norm) / fnorm(df2_noise)
X2_035_norm = df2norm + 0.35 * df2_noise * fnorm(df2norm) / fnorm(df2_noise)

# Prepare CMTF input
Z_010 = setupCMTFdata(list(X1_010@data, X2_010@data), modes, normalize=FALSE)
Z_025 = setupCMTFdata(list(X1_025@data, X2_025@data), modes, normalize=FALSE)
Z_035 = setupCMTFdata(list(X1_035@data, X2_035@data), modes, normalize=FALSE)

Z_010_norm = setupCMTFdata(list(X1_010_norm@data, X2_010_norm@data), modes, normalize=FALSE)
Z_025_norm = setupCMTFdata(list(X1_025_norm@data, X2_025_norm@data), modes, normalize=FALSE)
Z_035_norm = setupCMTFdata(list(X1_035_norm@data, X2_035_norm@data), modes, normalize=FALSE)

# Settings
nstart = 100
rel_tol = 1e-8
grad_tol = 1e-8
max_fn = 10^4
max_iter = 10^3
numCores = detectCores()
print(numCores)

CMTF_010 = cmtf_opt(Z_010, R, initialization="random", rel_tol=rel_tol, grad_tol=grad_tol, max_fn=max_fn, max_iter=max_iter, nstart=nstart, numCores=numCores, allOutput=TRUE)
CMTF_025 = cmtf_opt(Z_025, R, initialization="random", rel_tol=rel_tol, grad_tol=grad_tol, max_fn=max_fn, max_iter=max_iter, nstart=nstart, numCores=numCores, allOutput=TRUE)
CMTF_035 = cmtf_opt(Z_035, R, initialization="random", rel_tol=rel_tol, grad_tol=grad_tol, max_fn=max_fn, max_iter=max_iter, nstart=nstart, numCores=numCores, allOutput=TRUE)

CMTF_010_norm = cmtf_opt(Z_010_norm, R, initialization="random", rel_tol=rel_tol, grad_tol=grad_tol, max_fn=max_fn, max_iter=max_iter, nstart=nstart, numCores=numCores, allOutput=TRUE)
CMTF_025_norm = cmtf_opt(Z_025_norm, R, initialization="random", rel_tol=rel_tol, grad_tol=grad_tol, max_fn=max_fn, max_iter=max_iter, nstart=nstart, numCores=numCores, allOutput=TRUE)
CMTF_035_norm = cmtf_opt(Z_035_norm, R, initialization="random", rel_tol=rel_tol, grad_tol=grad_tol, max_fn=max_fn, max_iter=max_iter, nstart=nstart, numCores=numCores, allOutput=TRUE)

# Recreate input Fac using the component sizes to sort by size
inputFac = list(A, B, C, D, E)
inputFacnorm = list(Anorm, Bnorm, Cnorm, Dnorm, Enorm)

result_Tensor_Tensor = array(0L, c(6, length(CMTF_010)))
for(i in 1:length(CMTF_010)){
  result_Tensor_Tensor[1,i] = FMS_to_simulation(inputFac, CMTF_010[[i]]$Fac, modes)
  result_Tensor_Tensor[2,i] = FMS_to_simulation(inputFac, CMTF_025[[i]]$Fac, modes)
  result_Tensor_Tensor[3,i] = FMS_to_simulation(inputFac, CMTF_035[[i]]$Fac, modes)

  result_Tensor_Tensor[4,i] = FMS_to_simulation(inputFacnorm, CMTF_010_norm[[i]]$Fac, modes)
  result_Tensor_Tensor[5,i] = FMS_to_simulation(inputFacnorm, CMTF_025_norm[[i]]$Fac, modes)
  result_Tensor_Tensor[6,i] = FMS_to_simulation(inputFacnorm, CMTF_035_norm[[i]]$Fac, modes)
}

result_Tensor_Tensor
rowSums((result_Tensor_Tensor > 0.99^5) / nstart) * 100
rowMeans(result_Tensor_Tensor)

saveRDS(inputFac, "./noise_sim_tensor_tensor_inputFac.RDS")
saveRDS(inputFacnorm, "./noise_sim_tensor_tensor_inputFacnorm.RDS")

saveRDS(Z_010, "./noise_sim_tensor_tensor_Z_010.RDS")
saveRDS(Z_025, "./noise_sim_tensor_tensor_Z_025.RDS")
saveRDS(Z_035, "./noise_sim_tensor_tensor_Z_035.RDS")
saveRDS(Z_010_norm, "./noise_sim_tensor_tensor_Z_010_norm.RDS")
saveRDS(Z_025_norm, "./noise_sim_tensor_tensor_Z_025_norm.RDS")
saveRDS(Z_035_norm, "./noise_sim_tensor_tensor_Z_035_norm.RDS")

saveRDS(CMTF_010, "./noise_sim_tensor_tensor_CMTF_010.RDS")
saveRDS(CMTF_025, "./noise_sim_tensor_tensor_CMTF_025.RDS")
saveRDS(CMTF_035, "./noise_sim_tensor_tensor_CMTF_035.RDS")
saveRDS(CMTF_010_norm, "./noise_sim_tensor_tensor_CMTF_010_norm.RDS")
saveRDS(CMTF_025_norm, "./noise_sim_tensor_tensor_CMTF_025_norm.RDS")
saveRDS(CMTF_035_norm, "./noise_sim_tensor_tensor_CMTF_035_norm.RDS")

write.csv(result_Tensor_Tensor, "./noise_sim_tensor_tensor_output.csv")
