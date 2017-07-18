require(parallel)

## Vectorized upper triangle of the robust correlation matrix
cormad.vec <- function(data, jobs=15){
  p     <- ncol(data)
  idx   <- which(upper.tri(matrix(0L, ncol=p, nrow=p ), diag=FALSE), arr.ind=TRUE)
  ncors <- { p * {p-1} / 2 }
  ans   <- rep(0, times = ncors ) 
  
  ## Constants
  cost   <- 1.4826
  sqrt2  <- sqrt(2) 
  
  ## Compute principal variables 
  for (k in 1:p){
    m    <- median(data[,k] )
    s    <- cost * median( abs( data[,k] -   m ) )
    data[,k] <- { data[,k] - m} / { sqrt2 * s }
  }
  
  ## Compute pairwise correlations 
  ans <- unlist(mclapply(1:ncors, FUN=function(k) {
    i      <- idx[k,1]  ## row index of the upper triangle 
    j      <- idx[k,2]  ## col index of the upper triangle
    U      <- data[,i] + data[,j]
    V      <- data[,i] - data[,j]
    mad_U2 <- { cost * median( abs( U - median(U) )  ) }^2
    mad_V2 <- { cost * median( abs( V - median(V) )  ) }^2
    { mad_U2 - mad_V2 }  /  {mad_U2  + mad_V2} 
  }, mc.cores=jobs))
  
  return(ans)
}


##  cormad vector from <cormad.vec>  to cormad  matrix conversion 
cormad.vec2mat <- function(rho){
  
  ncors <- length(rho)
  p     <- {1+sqrt(1 + 4*2*ncors )}/2
  ans   <- matrix(1, ncol=p, nrow=p)
  idx   <- which(upper.tri(ans, diag=FALSE), arr.ind=TRUE)
  
  for (k in 1:ncors ){
    i <- idx[k,1]  ## row index of the upper triangle 
    j <- idx[k,2]  ## col index of the upper triangle
    ans[i,j]<-ans[j,i]<- rho[k]
  }
  
  return(ans)
}
