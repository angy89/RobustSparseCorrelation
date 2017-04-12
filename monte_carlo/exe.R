library(devtools)
source('R_functions.R')

## Global paramters
R         <- 250
nsplits   <- 25
p         <- 500
cont.rate <- 0.1

##  Correlation models
CorModel      <- array(0, dim=c(p,p,2))
## Spherical
CorModel[,,1] <- diag(1, p)
## Block diagonal
p1 <-  round(p/2)
tmp <- diag(1, p)
tmp[1:p1, 1:p1] <- cov.ar(p1, rho = 0.7, marg.vars=NULL)
CorModel[,,2] <- tmp

## Load setups
setups    <- read.table("setups.dat", header=TRUE, stringsAsFactors=FALSE)

## Do cases
for (ncase in 1:ncol(setups)) {
   set.seed(ncase * 19051977)
   source('docase.R')
}
