# Robust Sparse Correlation

This is the R implementation of a robust correlation matrix estimator usefull in case of high-dimensional data to remove the bias. 

It also implement an adaptive thresholding technique based on cross validation tecnique to regularize the correlation matrix.

### Demo
```R
source("cormad.R")
source("thresholding.R")

#load sample data
library("grndata")
data = syntren300.data
gold = syntren300.net

#evaluate robust correlation
R_vec = cormad.vec(red_data)
R = cormad.vec2mat(R_vec)

#compute thresholding statistics
FLOSSES = thresholding(data = red_data,nsplits = nsplits,n_cores = n_cores)
colSumsFlosses = rbind(colSumsFlosses,colSums(FLOSSES))

#find the threshold used to cut the correlation matrix
th = which.min(colSums(FLOSSES))/100

#cut the matrix
R = R[abs(R)<th]=0
```

### Simulation Code

The monte_carlo folder contains source code for simulation. The simulation can be executed with the following R command:

```R
source('exe.R')

