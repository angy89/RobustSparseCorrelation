## s_pred and s_true are two binary adjacency matrices. 
## s_pred is the prediceted matrix, while s_true is the gold standard
compute_statistics=function(s_pred,s_true){
  s_tp = sum(s_pred * s_true)
  s_tn = sum((1 - s_pred) * (1 - s_true))
  s_fp = sum(s_pred * (1 - s_true))
  s_fn = sum(s_true * (1 - s_pred))
  
  se = s_tp / (s_tp + s_fn)
  sp = s_tn / (s_tn + s_fp)
  acc = (s_tp  + s_tn)/ (s_tp + s_fp + s_fn + s_tn)
  f1_score = (2 * s_tp)/(2*s_tp + s_fp + s_fn)
  toRet = list('sensitivity'=se,'specificity'=sp,'accuracy'=acc,'f1score'=f1_score)
  return(toRet)
}

compare_statistics_after_threshold = function(listMat=NULL,gold = NULL){
  library(minet)
  
  if(is.null(gold)){
    print(" Error: Pleas provide the gold standard network! \n")
    return(-1)
  }
  
  if(is.null(listMat)){
    print(" Error: At least one between the computed network must be provided! \n")
    return(-1)
  }
  
  M = matrix(0,4,length(listMat))
  rownames(M)=c("Specificity","Sensitivity","Accuracy","F1")
  
  if(is.null(names(listMat))){
    print(" Warning: Correlation techniques not specified! \n")
  }else{
    colnames(M) = names(listMat)#c("RSC","TBW","Pearson","Spearman","Kendall")
  }
  
  for(i in 1:length(listMat)){
    MM = listMat[[i]]
    MM[is.na(MM)] = 0
    R_th_stats = compute_statistics(s_pred = aracne(mim = abs(MM)^2),s_true = gold)
    M["Specificity",i] = R_th_stats$specificity
    M["Sensitivity",i] = R_th_stats$sensitivity
    M["Accuracy",i] = R_th_stats$accuracy
    M["F1",i] = R_th_stats$f1score
  }

  return(M)
  
}
