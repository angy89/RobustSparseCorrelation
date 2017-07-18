## Load case parameters
##
ModelNum  <- setups[ncase, 2]
ModelName <- setups[ncase, 3]
cr        <- setups[ncase, 4]
p         <- floor(n * cr)
fn        <- paste(formatC(ncase, width=2, format = "d", flag = "0"))



## Build the correlation models
##
CorModel      <- array(0, dim=c(p,p,2))
## Spherical
CorModel[,,1] <- diag(1, p)
## Block diagonal
p1 <-  round(p/2)
tmp <- diag(1, p)
tmp[1:p1, 1:p1] <- cov.ar(p1, rho = 0.7, marg.vars=NULL)
CorModel[,,2] <- tmp



## Pick the correlation model
Sig      <- CorModel[ , , ModelNum]
eigenSig <- eigen(Sig, symmetric = TRUE)


## Initialize containers
Mets  <- c('Rn',  'Tau',  'Rho',  'Tbw',  'RMAD',
           'TRn', 'TTau', 'TRho', 'TTbw', 'RSC')
FROB  <- matrix(NA, nrow=R, ncol=10)
colnames(FROB) <- Mets
EIGL1 <- FROB
##
FITSCOUNT <- array(NA, dim=c(R,10,4)) ## FITSCOUNT[R, estimator, count_type]
dimnames(FITSCOUNT)[[2]] <- Mets
dimnames(FITSCOUNT)[[3]] <- c("zero2nonzero", "E1",    "E2",   "SS")
##
THRESHOLDS <- matrix(NA, nrow=R, ncol=5)
colnames(THRESHOLDS) <- Mets[6:10]




## Objects to be saved
Obj2Save <- c('FROB', 'EIGL1', 'FITSCOUNT', 'THRESHOLDS', 'last.rep')


cat("\n\n\n")
t0 <- Sys.time()


for (r in 1:R){
  message('Case #', ncase, ': ',  ModelName, ' model - ', 'Estimation at Replica ', r, ' of ', R)

  sim <- dgp(n=n, eigcov=eigenSig)

  ## Compute Rn
  nMet    <- 1L
  message(paste('........ Computing ', Mets[nMet], sep=''))
  X         <- Rn <- cor(sim$data)
  X_values  <- eigen(X, symmetric=TRUE)$values
  ## Store results
  FROB[r, nMet]      <- norm(X   - Sig, "F")
  EIGL1[r, nMet]     <- sum(abs(X_values - eigenSig$values))
  FITSCOUNT[r,nMet,] <- fitcounts(ref=Sig, fitted=X)


  ## Compute Tau
  nMet    <- 1L + nMet
  message(paste('........ Computing ', Mets[nMet], sep=''))
  X         <- Tau <- cor(sim$data, method = "kendall")
  X_values  <- eigen(X, symmetric=TRUE)$values
  ## Store results
  FROB[r, nMet]      <- norm(X   - Sig, "F")
  EIGL1[r, nMet]     <- sum(abs(X_values - eigenSig$values))
  FITSCOUNT[r,nMet,] <- fitcounts(ref=Sig, fitted=X)



  ## Compute Rho
  nMet    <- 1L + nMet
  message(paste('........ Computing ', Mets[nMet], sep=''))
  X         <- Rho <- cor(sim$data, method = 'spearman')
  X_values  <- eigen(X, symmetric=TRUE)$values
  ## Store results
  FROB[r, nMet]      <- norm(X   - Sig, "F")
  EIGL1[r, nMet]     <- sum(abs(X_values - eigenSig$values))
  FITSCOUNT[r,nMet,] <- fitcounts(ref=Sig, fitted=X)



  ## Compute TBW
  nMet    <- 1L + nMet
  ##
  ## @Angela: qui hai il TBW, tieni conto degli NA
  ##          se non vuoi eseguire il thresholding su TBW commenta il seguente
  ##          blocco di codice
  ##
  message(paste('........ Computing ', Mets[nMet], sep=''))
  vecTbw    <- TBW.vec(sim$data)
  X         <- Tbw  <- vec2mat(vecTbw)
  X_values  <- eigen(X, symmetric=TRUE)$values
  ## Store results
  FROB[r, nMet]      <- norm(X   - Sig, "F")
  EIGL1[r, nMet]     <- sum(abs(X_values - eigenSig$values))
  FITSCOUNT[r,nMet,] <- fitcounts(ref=Sig, fitted=X)



  ## Compute RMAD
  nMet    <- 1L + nMet
  message(paste('........ Computing ', Mets[nMet], sep=''))
  vecRmad   <- cormad.vec(sim$data)
  X         <- Rmad  <- vec2mat(vecRmad)
  X_values  <- eigen(X, symmetric=TRUE)$values
  ## Store results
  FROB[r, nMet]      <- norm(X   - Sig, "F")
  EIGL1[r, nMet]     <- sum(abs(X_values - eigenSig$values))
  FITSCOUNT[r,nMet,] <- fitcounts(ref=Sig, fitted=X)



  ## Compute TRn
  nMet    <- 1L + nMet
  message(paste('........ Computing ', Mets[nMet], sep=''))
  hsample <- THRESHOLDS[r, {nMet-5}] <-  cv.Rn(data=sim$data, nsplits=nsplits)$threshold
  X <- Rn
  X[abs(X) <= hsample] <- 0
  X_values  <- eigen(X, symmetric=TRUE)$values
  ## Store results
  FROB[r, nMet]      <- norm(X   - Sig, "F")
  EIGL1[r, nMet]     <- sum(abs(X_values - eigenSig$values))
  FITSCOUNT[r,nMet,] <- fitcounts(ref=Sig, fitted=X)



  ## Compute TTau
  nMet    <- 1L + nMet
  message(paste('........ Computing ', Mets[nMet], sep=''))
  hsample <- THRESHOLDS[r, {nMet-5}] <- cv.Tau(data=sim$data, nsplits=nsplits)$threshold
  X <- Tau
  X[abs(X) <= hsample] <- 0
  X_values  <- eigen(X, symmetric=TRUE)$values
  ## Store results
  FROB[r, nMet]      <- norm(X   - Sig, "F")
  EIGL1[r, nMet]     <- sum(abs(X_values - eigenSig$values))
  FITSCOUNT[r,nMet,] <- fitcounts(ref=Sig, fitted=X)




  ## Compute TRho
  nMet    <- 1L + nMet
  message(paste('........ Computing ', Mets[nMet], sep=''))
  hsample <- THRESHOLDS[r, {nMet-5}] <- cv.Rho(data=sim$data, nsplits=nsplits)$threshold
  X <- Rho
  X[abs(X) <= hsample] <- 0
  X_values  <- eigen(X, symmetric=TRUE)$values
  ## Store results
  FROB[r, nMet]      <- norm(X   - Sig, "F")
  EIGL1[r, nMet]     <- sum(abs(X_values - eigenSig$values))
  FITSCOUNT[r,nMet,] <- fitcounts(ref=Sig, fitted=X)





  ## Compute TTBW
  nMet    <- 1L + nMet
  ##
  ## @Angela: qui hai il TBW tagliato, tieni conto degli NA
  ##          se non vuoi eseguire il thresholding su TBW commenta il seguente
  ##          blocco di codice
  ##
  message(paste('........ Computing ', Mets[nMet], sep=''))
  hsample <- THRESHOLDS[r, {nMet-5}] <- cv.TBW(data=sim$data, nsplits=nsplits)$threshold
  X <- Tbw
  X[abs(X) <= hsample] <- 0
  X_values  <- eigen(X, symmetric=TRUE)$values
  ## Store results
  FROB[r, nMet]      <- norm(X   - Sig, "F")
  EIGL1[r, nMet]     <- sum(abs(X_values - eigenSig$values))
  FITSCOUNT[r,nMet,] <- fitcounts(ref=Sig, fitted=X)




  ## Compute RSC
  nMet    <- 1L + nMet
  message(paste('........ Computing ', Mets[nMet], sep=''))
  hsample <- THRESHOLDS[r, {nMet-5}] <- cv.Rmad(data=sim$data, nsplits=nsplits)$threshold
  X <- Rmad
  X[abs(X) <= hsample] <- 0
  X_values  <- eigen(X, symmetric=TRUE)$values
  ## Store results
  FROB[r, nMet]      <- norm(X   - Sig, "F")
  EIGL1[r, nMet]     <- sum(abs(X_values - eigenSig$values))
  FITSCOUNT[r,nMet,] <- fitcounts(ref=Sig, fitted=X)


   last.rep <- r


  if(r>=10 & {(r%%10)==0}){
    system(paste("rm   ",  fn, ".RData", sep=""))
    save(list=Obj2Save, file=paste(fn, ".RData", sep=""))
  }

  ## Compute remaining time
  dT <- difftime(Sys.time(), t0, units="secs")
  pretime <- as.numeric(dT/last.rep) * {R-last.rep}  ## predicted time to finish in secs
  cat("\n\n\n")
  message("Predicted time to finish:.....", round(pretime/{3600},5), " (hours)")
  message("Predicted end date:...........", Sys.time()+pretime)
  cat("\n\n\n")
}



## END
cat("\n\n")
dT <- difftime(Sys.time(), t0, units="secs")
message("Finished:.....................", Sys.time())
message("Total time taken (hours):.....", round(dT/3600,5))
cat("\n\n")


SessionInfo <- session_info()
EndDate     <- Sys.time()
Obj2Save    <- c(Obj2Save, 'SessionInfo', 'EndDate')
save(list=Obj2Save, file=paste(fn, ".RData", sep=""))


