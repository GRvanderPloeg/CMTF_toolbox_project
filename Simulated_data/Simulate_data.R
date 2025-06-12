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

# Save loadings
allLoadings = list(scores, loadings1, timeLoadings1, loadings2, timeLoadings2, loadings3, timeLoadings3)
saveRDS(allLoadings, "./Simulated_data/input_loadings.RDS")
