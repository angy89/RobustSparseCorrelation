################################################################################
## Software for cross-validation of Rn, Rho, Tau
################################################################################


## Cross-validation for sample correlation
cv.Rn <- function(data, nsplits=10, monitor=1, isPar = TRUE, n_cores = 8){
  n <- nrow(data)
  p <- ncol(data)
  tgrid    <- seq(0.005, 1, by=0.01)
  ngrid    <- length(tgrid)
  n1       <- n - floor(n/log(n))
  n2       <- n - n1
  
  ## Estimate expected loss
  if(isPar){
    registerDoMC(n_cores)
    FLOSSES <- foreach(i=1:nsplits, .combine='rbind') %dopar% {
      if (monitor>=1) { message("Split: ", i, " (of ", nsplits, ")") }
      idx <- sample(1:n, size=n1, replace=FALSE)
      C1 <- cor(  data[ idx,  ] )
      C2 <- cor(  data[-idx,  ] )
      
      foreach(k=1:ngrid, .combine='c') %do% {
        C1[  abs(C1) <  tgrid[k]  ] <- 0
        diag(C1) <- rep(1,p)
        norm(C1-C2, "F")
      }
    }  
  }else{
    #C1 <- C2 <- matrix(0, nrow=p, ncol=p)
    FLOSSES  <- matrix(0, nrow=nsplits, ncol=ngrid)
    
    for (i in 1:nsplits){
      if (monitor>=1) { message("Split: ", i, " (of ", nsplits, ")") }
      idx <- sample(1:n, size=n1, replace=FALSE)
      C1 <- cor(  data[ idx,  ]  )
      C2 <- cor(  data[-idx,  ]  )
      for (k in 1:ngrid){
        C1[  abs(C1) <  tgrid[k]  ] <- 0
        diag(C1) <- rep(1,p)
        FLOSSES[i,k] <- norm(C1-C2, "F")
      }
    }
  }
  
  avgLoss <- colMeans(FLOSSES)
  hstar   <- tgrid[which.min(avgLoss)]
  
  ans <- list()
  ans$avgloss   <- avgLoss
  ans$threshold <- hstar
  
  return(ans)
}

## Cross-validation for Kendall Tau
cv.Tau <- function(data, nsplits=10, monitor=1,isPar = TRUE, n_cores = 32){
  n <- nrow(data)
  p <- ncol(data)
  tgrid    <- seq(0.005, 1, by=0.01)
  ngrid    <- length(tgrid)
  n1       <- n - floor(n/log(n))
  n2       <- n - n1 
  
  ## Estimate expected loss
  if(isPar){
    registerDoMC(n_cores)
    FLOSSES <- foreach(i=1:nsplits, .combine='rbind') %dopar% {
      if (monitor>=1) { message("Split: ", i, " (of ", nsplits, ")") }
      idx <- sample(1:n, size=n1, replace=FALSE)
      C1 <- cor(  data[ idx,  ], method = 'kendall')
      C2 <- cor(  data[-idx,  ], method = 'kendall')
      
      foreach(k=1:ngrid, .combine='c') %do% {
        C1[  abs(C1) <  tgrid[k]  ] <- 0
        diag(C1) <- rep(1,p)
        norm(C1-C2, "F")
      }
    }  
  }else{
    #C1 <- C2 <- matrix(0, nrow=p, ncol=p)
    FLOSSES  <- matrix(0, nrow=nsplits, ncol=ngrid)
    
    for (i in 1:nsplits){
      if (monitor>=1) { message("Split: ", i, " (of ", nsplits, ")") }
      idx <- sample(1:n, size=n1, replace=FALSE)
      C1 <- cor(  data[ idx,  ], method = 'kendall'  )
      C2 <- cor(  data[-idx,  ], method = 'kendall'  )
      for (k in 1:ngrid){
        C1[  abs(C1) <  tgrid[k]  ] <- 0
        diag(C1) <- rep(1,p)
        FLOSSES[i,k] <- norm(C1-C2, "F")
      }
    }
  }
  
  avgLoss <- colMeans(FLOSSES)
  hstar   <- tgrid[which.min(avgLoss)]
  
  ans <- list()
  ans$avgloss   <- avgLoss
  ans$threshold <- hstar
  
  return(ans)
}

## Cross-validation for Spearman Rho
cv.Rho <- function(data, nsplits=10, monitor=1,isPar = TRUE, n_cores = 8){
  n <- nrow(data)
  p <- ncol(data)
  tgrid    <- seq(0.005, 1, by=0.01)
  ngrid    <- length(tgrid)
  n1       <- n - floor(n/log(n))
  n2       <- n - n1
  
  ## Estimate expected loss
  if(isPar){
    registerDoMC(n_cores)
    FLOSSES <- foreach(i=1:nsplits, .combine='rbind') %dopar% {
      if (monitor>=1) { message("Split: ", i, " (of ", nsplits, ")") }
      idx <- sample(1:n, size=n1, replace=FALSE)
      C1 <- cor(  data[ idx,  ], method = 'spearman')
      C2 <- cor(  data[-idx,  ], method = 'spearman')
      
      foreach(k=1:ngrid, .combine='c') %do% {
        C1[  abs(C1) <  tgrid[k]  ] <- 0
        diag(C1) <- rep(1,p)
        norm(C1-C2, "F")
      }
    }  
  }else{
    #C1 <- C2 <- matrix(0, nrow=p, ncol=p)
    FLOSSES  <- matrix(0, nrow=nsplits, ncol=ngrid)
    
    for (i in 1:nsplits){
      if (monitor>=1) { message("Split: ", i, " (of ", nsplits, ")") }
      idx <- sample(1:n, size=n1, replace=FALSE)
      C1 <- cor(  data[ idx,  ], method = 'spearman'  )
      C2 <- cor(  data[-idx,  ], method = 'spearman'  )
      for (k in 1:ngrid){
        C1[  abs(C1) <  tgrid[k]  ] <- 0
        diag(C1) <- rep(1,p)
        FLOSSES[i,k] <- norm(C1-C2, "F")
      }
    }
  }
  avgLoss <- colMeans(FLOSSES)
  hstar   <- tgrid[which.min(avgLoss)]
  
  ans <- list()
  ans$avgloss   <- avgLoss
  ans$threshold <- hstar
  
  return(ans)
}

