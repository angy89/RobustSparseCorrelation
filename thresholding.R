#Data driven thresholding for correlation matrix

thresholding = function(data,nsplits = 1000,tgrid = seq(0.005, 1, by=0.01),n_cores=1,monitor=1){
  ## Cross-validation    
  n        <- nrow(data)    ## patients
  p        <- ncol(data)    ## genes
  nc       <- p * {p-1} / 2
  nsplits  <- nsplits
  tgrid    <- tgrid
  ngrid    <- length(tgrid)
  n1       <- n - floor(n/log(n))
  n2       <- n-n1
  C1 <- C2 <- rep(0, times=nc) 
  FLOSSES  <- matrix(0, nrow=nsplits, ncol=ngrid)
  
  
  ## Estimate expected loss    
  for (i in 1:nsplits){
    if (monitor>=1) { message("Split: ", i, " (of ", nsplits, ")") }
   # message("Split: ", i, " (of ", nsplits, ")")
    idx <- sample(1:n, size=n1, replace=FALSE)
    C1 <- cormad.vec(  data[ idx,  ],n_cores  ) 
    C2 <- cormad.vec(  data[-idx,  ],n_cores  )
    
    if(sum(is.na(C1))>0 || sum(is.na(C2))>0)
      print("NA in RSC...\n")
    
    for (k in 1:ngrid){
      C1[  abs(C1) <  tgrid[k]  ] <- 0
      FLOSSES[i,k] <- sum( 2 * { C1 - C2}^2 ) 
    }
    
  }
  
  avgLoss <- colMeans(FLOSSES)
  #hstar   <- tgrid[which.min(avgLoss)]
  
  ans = list()
  ans$avgloss   <- avgLoss
  ans$threshold <- which.min(colSums(FLOSSES))/100
  
  
  return(ans)
}
