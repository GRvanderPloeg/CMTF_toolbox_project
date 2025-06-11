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

  # Calculate FMS
  FMS_result = 1:numComponents
  for(i in 1:numComponents){
    proposedFMS = 1:numComponents
    for(j in 1:numComponents){
      FMS = 1
      ksi = 0
      ksi_hat = 0

      for(k in 1:numModes){
        vect1 = normalizedFac$Fac[[k]][,i]
        vect2 = normalizedInputFac$Fac[[k]][,j]
        FMS = FMS * ((t(vect2) %*% vect1))
      }

      for(p in 1:numDatasets){
        ksi_hat = ksi_hat + normalizedFac$normsPerDataset[p,i]
        ksi = ksi + normalizedInputFac$normsPerDataset[p,j]
      }
      ksi_term = ((abs(ksi - ksi_hat)) / max(ksi, ksi_hat))
      proposedFMS[j] = (1 - ksi_term) * abs(FMS)
    }
    FMS_result[i] = max(proposedFMS)
  }
  return(min(FMS_result))
}

R = 3
I = 75
J = 75
K = 75
L = 75
M = 75
modes = list(c(1,2), c(1,3,4), c(4,5))

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
df1 = as.tensor(reinflateMatrix(A, B))
df2 = as.tensor(reinflateTensor(A, C, D))
df3 = as.tensor(reinflateMatrix(D, E))

df1norm = as.tensor(reinflateMatrix(Anorm, Bnorm))
df2norm = as.tensor(reinflateTensor(Anorm, Cnorm, Dnorm))
df3norm = as.tensor(reinflateMatrix(Dnorm, Enorm))

# Generate noise
df1_noise = as.tensor(array(rnorm(I*J), c(I,J)))
df2_noise = as.tensor(array(rnorm(I*K*L), c(I,K,L)))
df3_noise = as.tensor(array(rnorm(L*M), c(L,M)))

# Create final datasets by adding the noise of different sizes to the generated data
X1_010 = df1 + 0.10 * df1_noise * fnorm(df1) / fnorm(df1_noise)
X1_025 = df1 + 0.25 * df1_noise * fnorm(df1) / fnorm(df1_noise)
X1_035 = df1 + 0.35 * df1_noise * fnorm(df1) / fnorm(df1_noise)

X2_010 = df2 + 0.10 * df2_noise * fnorm(df2) / fnorm(df2_noise)
X2_025 = df2 + 0.25 * df2_noise * fnorm(df2) / fnorm(df2_noise)
X2_035 = df2 + 0.35 * df2_noise * fnorm(df2) / fnorm(df2_noise)

X3_010 = df3 + 0.10 * df3_noise * fnorm(df3) / fnorm(df3_noise)
X3_025 = df3 + 0.25 * df3_noise * fnorm(df3) / fnorm(df3_noise)
X3_035 = df3 + 0.35 * df3_noise * fnorm(df3) / fnorm(df3_noise)

X1_010_norm = df1norm + 0.10 * df1_noise * fnorm(df1norm) / fnorm(df1_noise)
X1_025_norm = df1norm + 0.25 * df1_noise * fnorm(df1norm) / fnorm(df1_noise)
X1_035_norm = df1norm + 0.35 * df1_noise * fnorm(df1norm) / fnorm(df1_noise)

X2_010_norm = df2norm + 0.10 * df2_noise * fnorm(df2norm) / fnorm(df2_noise)
X2_025_norm = df2norm + 0.25 * df2_noise * fnorm(df2norm) / fnorm(df2_noise)
X2_035_norm = df2norm + 0.35 * df2_noise * fnorm(df2norm) / fnorm(df2_noise)

X3_010_norm = df3norm + 0.10 * df3_noise * fnorm(df3norm) / fnorm(df3_noise)
X3_025_norm = df3norm + 0.25 * df3_noise * fnorm(df3norm) / fnorm(df3_noise)
X3_035_norm = df3norm + 0.35 * df3_noise * fnorm(df3norm) / fnorm(df3_noise)

# Prepare CMTF input
Z_010 = setupCMTFdata(list(X1_010@data, X2_010@data, X3_010@data), modes, normalize=FALSE)
Z_025 = setupCMTFdata(list(X1_025@data, X2_025@data, X3_025@data), modes, normalize=FALSE)
Z_035 = setupCMTFdata(list(X1_035@data, X2_035@data, X3_035@data), modes, normalize=FALSE)

Z_010_norm = setupCMTFdata(list(X1_010_norm@data, X2_010_norm@data, X3_010_norm@data), modes, normalize=FALSE)
Z_025_norm = setupCMTFdata(list(X1_025_norm@data, X2_025_norm@data, X3_025_norm@data), modes, normalize=FALSE)
Z_035_norm = setupCMTFdata(list(X1_035_norm@data, X2_035_norm@data, X3_035_norm@data), modes, normalize=FALSE)

# Settings
nstart = 100 # 30
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

result_Matrix_Tensor_Matrix = array(0L, c(6, length(CMTF_010)))
for(i in 1:length(CMTF_010)){
  result_Matrix_Tensor_Matrix[1,i] = FMS_to_simulation(inputFac, CMTF_010[[i]]$Fac, modes)
  result_Matrix_Tensor_Matrix[2,i] = FMS_to_simulation(inputFac, CMTF_025[[i]]$Fac, modes)
  result_Matrix_Tensor_Matrix[3,i] = FMS_to_simulation(inputFac, CMTF_035[[i]]$Fac, modes)

  result_Matrix_Tensor_Matrix[4,i] = FMS_to_simulation(inputFacnorm, CMTF_010_norm[[i]]$Fac, modes)
  result_Matrix_Tensor_Matrix[5,i] = FMS_to_simulation(inputFacnorm, CMTF_025_norm[[i]]$Fac, modes)
  result_Matrix_Tensor_Matrix[6,i] = FMS_to_simulation(inputFacnorm, CMTF_035_norm[[i]]$Fac, modes)
}

result_Matrix_Tensor_Matrix
rowSums((result_Matrix_Tensor_Matrix > 0.99^5) / nstart) * 100
rowMeans(result_Matrix_Tensor_Matrix)

saveRDS(inputFac, "./noise_sim_matrix_tensor_matrix_inputFac.RDS")
saveRDS(inputFacnorm, "./noise_sim_matrix_tensor_matrix_inputFacnorm.RDS")

saveRDS(Z_010, "./noise_sim_matrix_tensor_matrix_Z_010.RDS")
saveRDS(Z_025, "./noise_sim_matrix_tensor_matrix_Z_025.RDS")
saveRDS(Z_035, "./noise_sim_matrix_tensor_matrix_Z_035.RDS")
saveRDS(Z_010_norm, "./noise_sim_matrix_tensor_matrix_Z_010_norm.RDS")
saveRDS(Z_025_norm, "./noise_sim_matrix_tensor_matrix_Z_025_norm.RDS")
saveRDS(Z_035_norm, "./noise_sim_matrix_tensor_matrix_Z_035_norm.RDS")

saveRDS(CMTF_010, "./noise_sim_matrix_tensor_matrix_CMTF_010.RDS")
saveRDS(CMTF_025, "./noise_sim_matrix_tensor_matrix_CMTF_025.RDS")
saveRDS(CMTF_035, "./noise_sim_matrix_tensor_matrix_CMTF_035.RDS")
saveRDS(CMTF_010_norm, "./noise_sim_matrix_tensor_matrix_CMTF_010_norm.RDS")
saveRDS(CMTF_025_norm, "./noise_sim_matrix_tensor_matrix_CMTF_025_norm.RDS")
saveRDS(CMTF_035_norm, "./noise_sim_matrix_tensor_matrix_CMTF_035_norm.RDS")

write.csv(result_Matrix_Tensor_Matrix, "./noise_sim_matrix_tensor_matrix_output.csv")
