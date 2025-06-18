library(rTensor)
library(CMTFtoolbox)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(parallel)
library(doParallel)
library(foreach)

set.seed(123)

# Settings
numComponents = 8
cvFolds = 5
nstart = 10
normalize = TRUE
normY = 1
numCores = parallel::detectCores()
pis = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

# Load data
X1_final = readRDS("./Sim1_X1.RDS")
X2_final = readRDS("./Sim1_X2.RDS")
X3_final = readRDS("./Sim1_X3.RDS")
Y_final = as.matrix(readRDS("./Sim1_Y.RDS"))

# Process data
datasets = list(X1_final, X2_final, X3_final)
modes = list(c(1,2,3), c(1,4,5), c(1,6,7))
Z = CMTFtoolbox::setupCMTFdata(datasets, modes, normalize=TRUE)

numDatasets = length(datasets)
numModes = max(unlist(modes))
numSubjects = dim(datasets[[1]])[1]

# Create CV folds
indices = seq_len(numSubjects)
foldsPartition = split(indices, cut(seq_along(indices), breaks = cvFolds, labels = FALSE))
uniqueFolds = seq_len(cvFolds)

for(pi in pis){

  print(pi)

  # Save Ypred across the folds
  Ypred = matrix(NA, nrow=nrow(Y_final), ncol=nstart)

  for(i in uniqueFolds){

    print(i)
    testIdx = foldsPartition[[i]]
    trainIdx = setdiff(seq_len(numSubjects), testIdx)

    ## Prepare X
    Xtrain_final = list()
    Xtest_final = list()
    for(p in 1:numDatasets){

      if(length(dim(Z$object[[p]]@data))==3){
        Xtrain = rTensor::as.tensor(Z$object[[p]]@data[trainIdx, ,])
        Xtest = rTensor::as.tensor(Z$object[[p]]@data[testIdx, ,])
      } else{
        Xtrain = rTensor::as.tensor(Z$object[[p]]@data[trainIdx,])
        Xtest = rTensor::as.tensor(Z$object[[p]]@data[testIdx,])
      }

      # # Centering Xtrain
      # unfoldedXtrain = rTensor::k_unfold(Xtrain, 1)@data
      # means = colMeans(unfoldedXtrain, na.rm=TRUE)
      # unfoldedXtrain_cnt = sweep(unfoldedXtrain, 2, means, FUN="-")
      # Xtrain_cnt = rTensor::k_fold(unfoldedXtrain_cnt, m=1, modes=Xtrain@modes)
      #
      # # Use the means to center Xtest as well
      # unfoldedXtest = rTensor::k_unfold(Xtest, 1)@data
      # unfoldedXtest_cnt = sweep(unfoldedXtest, 2, means, FUN="-")
      # Xtest_cnt = rTensor::k_fold(unfoldedXtest_cnt, m=1, modes=Xtest@modes)
      #
      # # Scaling Xtrain
      # unfoldedXtrain = rTensor::k_unfold(Xtrain_cnt, 2)@data
      # stds = apply(unfoldedXtrain, 1, function(x){stats::sd(x, na.rm=TRUE)})
      # unfoldedXtrain_scl = sweep(unfoldedXtrain, 1, stds, FUN="/")
      # Xtrain_cnt_scl = rTensor::k_fold(unfoldedXtrain_scl, m=2, modes=Xtrain@modes)
      #
      # # Use the stds to scale Xtest as well
      # unfoldedXtest = rTensor::k_unfold(Xtest_cnt, 2)@data
      # unfoldedXtest_scl = sweep(unfoldedXtest, 1, stds, FUN="/")
      # Xtest_cnt_scl = rTensor::k_fold(unfoldedXtest_scl, m=2, modes=Xtest@modes)

      Xtrain_final = Xtrain
      Xtest_final = Xtest

    #   if(normalize){
    #     norm = rTensor::fnorm(Xtrain_cnt_scl)
    #     Xtrain_final[[p]] = Xtrain_cnt_scl@data / norm
    #     Xtest_final[[p]] = Xtest_cnt_scl@data / norm
    #   } else{
    #     Xtrain_final[[p]] = Xtrain_cnt_scl@data
    #     Xtest_final[[p]] = Xtest_cnt_scl@data
    #   }
    # }
    }

    Ztrain = setupCMTFdata(Xtrain_final, Z$modes, normalize=FALSE) # do not normalize again

    ## Prepare Y train
    Ytrain = Y_final[trainIdx]

    # Ymean = mean(Ytrain)
    # Ytrain_cnt = Ytrain - Ymean
    #
    # Ynorm = norm(Ytrain_cnt, "2")
    # Ytrain_normalized = Ytrain_cnt / Ynorm
    #
    # Ytrain_final = Ytrain_normalized * normY
    Ytrain_final = Ytrain

    # Each replicate is fitted with a single initialization and one core.
    models = CMTFtoolbox::acmtfr_opt(Ztrain,
                                    Ytrain_final,
                                    numComponents = 8,
                                    pi            = pi,
                                    method        = "L-BFGS",
                                    nstart        = nstart,
                                    numCores      = numCores,
                                    allOutput     = TRUE)

    # Predict Y for the test set data
    for(j in 1:nstart){
      model = models[[j]]
      pred = npred(model, list(Xtest_final[[1]], Xtest_final[[2]], Xtest_final[[3]]), Ztrain)
      # pred_norm = pred / normY
      # pred_cnt = pred_norm * Ynorm
      # pred_original = pred_cnt + Ymean
      pred_original = pred

      Ypred[testIdx,j] = pred_original
    }
  }

  saveRDS(Ypred, paste0("ACMTFR_CV_pi_", pi, ".RDS"))
}
