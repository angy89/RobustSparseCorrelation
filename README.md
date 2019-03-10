# Robust Sparse Correlation

This is the R implementation of a robust correlation matrix estimator usefull in case of high-dimensional data to remove the bias. 

It also implement an adaptive thresholding technique based on cross validation tecnique to regularize the correlation matrix.

Reference Paper:

> Serra, A., Coretto, P., Fratello, M., & Tagliaferri, R. (2017). Robust and sparse correlation matrix estimation for the analysis of high-dimensional genomics data. Bioinformatics, 34(4), 625-634.


### Demo
```R
source("cormad.R")
source("thresholding.R")

#load sample data
library("grndata")
data = syntren300.data
gold = syntren300.net
sample_size = 60
n_cores = 5

chips = sample(x = 1:nrow(data),size = sample_size,replace = FALSE)
red_data = data[chips,]

print("Evaluating RSC Correlation...\n")
RSC = thresholding(data = red_data,nsplits = nsplits,n_cores = n_cores,monitor = 0)

R_vec = cormad.vec(red_data,jobs = n_cores)
R = cormad.vec2mat(R_vec)
R_bin = bin_matrix(R,RSC$threshold)
```

### Simulation Code

The monte_carlo folder contains source code for simulation. The simulation can be executed with the following R command:

```R
source('exe.R')

