ModelNum  <- setups[ncase, 2]
ModelName <- setups[ncase, 3]
cr        <- setups[ncase, 4]
n         <- ceiling(p/cr)

## File naming
fn <- paste(formatC(ncase, width=2, format = "d", flag = "0"))

## Pick the correlation model
Sig      <- CorModel[ , , ModelNum]
eigenSig <- eigen(Sig, symmetric = TRUE)

## Initialize containers
FROB <- matrix(0, nrow=R, ncol=4)
colnames(FROB) <- c('SampleCor', 'ThresholdSampleCor', 'RobustCor', 'ThresholdRobustCor')
EIGL1 <- FROB
FITSCOUNT <- array(0, dim=c(R,4,4)) ## FITSCOUNT[R, estimator, count_type]
dimnames(FITSCOUNT)[[2]] <- as.list(colnames(FROB))
dimnames(FITSCOUNT)[[3]] <- c("zero2nonzero", "zero2high",    "high2zero",   "signswitch")


## Objects to be saved
Obj2Save <- c('FROB', 'EIGL1', 'FITSCOUNT', 'last.rep')


cat("\n\n\n")
t0 <- Sys.time()


for (r in 1:R){
  message('Case #', ncase, ': ',  ModelName, ' model - ', 'Estimation at Replica ', r, ' of ', R)

  sim <- dgp(n=n, eigcov=eigenSig, cont.rate=cont.rate)


  ## Sampling correlation
  message('........ Computing sample correlation')
  Rn        <- cor(sim$data)
  RnValues  <- eigen(Rn, symmetric=TRUE)$values
  ## Store results
  FROB[r, 1]      <- norm(Rn   - Sig, "F")
  EIGL1[r, 1]     <- sum(abs(RnValues - eigenSig$values))
  FITSCOUNT[r,1,] <- fitcounts(ref=Sig, fitted=Rn)




  ## Thresholded sample correlation
  message('........ Computing thresholded sample correlation')
  hsample <- cv.Rn(data=sim$data, nsplits=nsplits)$threshold
  TRn     <- Rn
  TRn[abs(Rn) <= hsample] <- 0
  TRnValues  <- eigen(TRn, symmetric=TRUE)$values
  ## Store results
  FROB[r, 2]      <- norm(TRn   - Sig, "F")
  EIGL1[r, 2]     <- sum(abs(TRnValues - eigenSig$values))
  FITSCOUNT[r,2,] <- fitcounts(ref=Sig, fitted=TRn)




  ## Robust correlation based on cormad
  message('........ Computing robust correlation')
  vecRmad <- cormad.vec(sim$data)
  Rmad    <- cormad.vec2mat(vecRmad)
  RmadValues <- eigen(Rmad, symmetric=TRUE)$values
  ## Store results
  FROB[r, 3]      <- norm(Rmad   - Sig, "F")
  EIGL1[r, 3]     <- sum(abs(RmadValues - eigenSig$values))
  FITSCOUNT[r,3,] <- fitcounts(ref=Sig, fitted=Rmad)



  ## RCS: robust sparse correlations
  message('........ Computing thresholded robust correlation')
  hmad <- cv.Rmad(data=sim$data, nsplits=nsplits)$threshold
  Rsc  <- Rmad
  Rsc[abs(Rmad) <= hmad] <- 0
  RscValues <- eigen(Rsc, symmetric=TRUE)$values
  ## Store results
  FROB[r, 4]      <- norm(Rsc   - Sig, "F")
  EIGL1[r, 4]     <- sum(abs(RscValues - eigenSig$values))
  FITSCOUNT[r,4,] <- fitcounts(ref=Sig, fitted=Rsc)

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









