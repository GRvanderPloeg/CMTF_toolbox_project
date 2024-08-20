library(CMTFtoolbox)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rTensor)
library(parallel)
library(multiway)

set.seed(123)

R = 3
I = 100
J = 100
K = 100
L = 100
modes = list(c(1,2,3), c(1,4))

# Generate loading matrices, normalize each vector to 1
A = array(rnorm(I*R), c(I, R))  # shared subject mode
B = array(rnorm(J*R), c(J, R))  # distinct feature mode of X1
C = array(rnorm(K*R), c(K, R))  # distinct condition mode of X1
D = array(rnorm(L*R), c(L, R))  # distinct feature mode of X2

# Generate datasets without noise
df1 = as.tensor(reinflateTensor(A, B, C))
df2 = as.tensor(reinflateMatrix(A, D))

# Set M% of data missing in df1 - missing completely at random (MCAR)
missingNess = c(seq(0.3, 0.9, 0.1), 0.95, 0.99)
missingIndices = list()
df1_list = list()
Zs = list()
for(i in 1:length(missingNess)){
  W = array(as.numeric(!array(runif(I*J*K), c(I,J,K)) < missingNess[i]), c(I,J,K))
  df1_missing = (W * df1)@data
  df1_missing[df1_missing == 0] = NA
  df1_list[[i]] = df1_missing
  missingIndices[[i]] = W
  Zs[[i]] = setupCMTFdata(list(df1_missing, df2@data), modes, normalize=FALSE)
}

# Create CP and CMTF models
cp_models = list()
cmtf_models = list()

for(i in 1:length(missingNess)){
  cp_models[[i]] = multiway::parafac(Zs[[i]]$object[[1]]@data, R, verbose=FALSE, nstart=100)
  cmtf_models[[i]] = cmtf_opt(Zs[[i]], R, initialization="random", nstart=100, numCores=parallel::detectCores())
}

# Calculate Matrix Completion diagnostic
TCS_CMTF = rep(0, length(missingNess))
TCS_CP = rep(0, length(missingNess))

for(i in 1:length(missingNess)){
  W = missingIndices[[i]]
  Xhat_CMTF = reinflateTensor(cmtf_models[[i]]$Fac[[1]], cmtf_models[[i]]$Fac[[2]], cmtf_models[[i]]$Fac[[3]])
  Xhat_CP = reinflateTensor(cp_models[[i]]$A, cp_models[[i]]$B, cp_models[[i]]$C)
  TCS_CMTF[i] = fnorm((1 - W) * (df1 - Xhat_CMTF)) / fnorm((1 - W) * df1)
  TCS_CP[i] = fnorm((1 - W) * (df1 - Xhat_CP)) / fnorm((1 - W) * df1)
}

# Save results
inputFac = list(A,B,C,D)
result = cbind(missingNess, TCS_CMTF, TCS_CP)

saveRDS(inputFac, "missing_sim_tensor_matrix_inputFac.RDS")
saveRDS(df1, "missing_sim_tensor_matrix_df1.RDS")
saveRDS(Zs, "missing_sim_tensor_matrix_Zs.RDS")
saveRDS(missingIndices, "missing_sim_tensor_matrix_Ws.RDS")
saveRDS(cmtf_models, "missing_sim_tensor_matrix_cmtf_models.RDS")
saveRDS(cp_models, "missing_sim_tensor_matrix_cp_models.RDS")
write.csv(result, "missing_sim_tensor_matrix_output.csv")
