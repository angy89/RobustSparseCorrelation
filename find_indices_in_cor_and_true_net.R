find_indices_in_cor_and_true_net = function(R,gold,th,main_text){
  
  gold_not_symm = gold
  gold = 0 * gold
  for(i in 1:nrow(gold_not_symm)){
    for(j in 1:ncol(gold_not_symm)){
      if(gold_not_symm[i,j]==1){
        gold[i,j] = gold[j,i] = 1
      }
    }  
  }
  
  true = gold[upper.tri(gold,diag=FALSE)]
  rcor = abs(R[upper.tri(R,diag = FALSE)])
  idx = order(rcor,decreasing = TRUE)
  rcor = rcor[idx]
  true = true[idx]
  true_pos = which(true==1)
  true_neg = which(true==0)
  
  plot(1,col="white",ylim = c(0,1),ylab="Correlation",xlab="",
       xlim = c(0,length(rcor)),main=main_text)
  #ßinvisible(sapply(X = true_neg,FUN = function(i)abline(v=i,col="green")))
  invisible(sapply(X = true_pos,FUN = function(i)abline(v=i,col="grey")))
  
  lines(rcor,col="black",lwd = 2)
  abline(h=th,col="black")

  return(list(rcor,true,true_pos))
  
}