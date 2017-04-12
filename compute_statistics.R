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

# evaluate sensitivity and specificity beteween the gold standard and the networks
# generated with aracne using the robust, the pearson, spearman and kendall correlations
compare_statistics = function(R_bin,pears,spearman,kendall,gold){
  library(minet)
  R_stats = compute_statistics(s_pred = aracne(mim = abs(R_bin)^2),s_true = gold)
  
  rob_spec = c()
  pears_spec = c()
  spearman_spec = c()
  kendall_spec = c()
  
  rob_sens = c()
  pears_sens = c()
  spearman_sens = c()
  kendall_sens = c()
  
  rob_acc = c()
  pears_acc = c()
  spearman_acc = c()
  kendall_acc = c()
  
  rob_f1 = c()
  pears_f1 = c()
  spearman_f1 = c()
  kendall_f1 = c()
  
  for(thi in seq(0,0.9,0.1)){
    pears_bin = pears
    pears_bin[abs(pears)<=thi] = 0
    pears_bin[abs(pears)>thi] = 1
    
    spearman_bin = spearman
    spearman_bin[abs(spearman)<=thi] = 0
    spearman_bin[abs(spearman)>thi] = 1
    
    kendall_bin = kendall
    kendall_bin[abs(kendall)<=thi] = 0
    kendall_bin[abs(kendall)>thi] = 1
    
    R_bin_th = R
    R_bin_th[abs(R)<=thi] = 0
    R_bin_th[abs(R)>thi] = 1
    
    R_th_stats = compute_statistics(s_pred = aracne(mim = abs(R_bin_th)),s_true = gold)
    pears_stats = compute_statistics(s_pred = aracne(mim = abs(pears_bin)^2),s_true = gold)
    spearman_stats = compute_statistics(s_pred = aracne(mim = abs(spearman_bin)^2),s_true = gold)
    kendall_stats = compute_statistics(s_pred = aracne(mim = abs(kendall_bin)^2),s_true = gold)
    
    rob_spec = c(rob_spec,R_th_stats$specificity)
    pears_spec = c(pears_spec,pears_stats$specificity)
    spearman_spec = c(spearman_spec,spearman_stats$specificity)
    kendall_spec = c(kendall_spec,kendall_stats$specificity)
    
    rob_sens = c(rob_sens,R_th_stats$sensitivity)
    pears_sens = c(pears_sens,pears_stats$sensitivity)
    spearman_sens= c(spearman_sens,spearman_stats$sensitivity)
    kendall_sens = c(kendall_sens,kendall_stats$sensitivity)
    
    rob_acc = c(rob_acc,R_th_stats$accuracy)
    pears_acc= c(pears_acc,pears_stats$accuracy)
    spearman_acc = c(spearman_acc,spearman_stats$accuracy)
    kendall_acc = c(kendall_acc,kendall_stats$accuracy)
    
    rob_f1 = c(rob_f1,R_th_stats$f1score)
    pears_f1= c(pears_f1,pears_stats$f1score)
    spearman_f1 = c(spearman_f1,spearman_stats$f1score)
    kendall_f1 = c(kendall_f1,kendall_stats$f1score)
  }
  
  line_type = c(1,3,5,6)
  seq_ = seq(0,0.9,0.1)
  lwd_ =2
  colors = c(1,3,5,6)
  
  plot(y=rep(1,length(seq_)),lwd = lwd_,x =seq_, ylab = "", col = "white",xaxt="n",
       xlab = "Thresholds",lty = line_type[1],ylim=c(0,1),type="l")
  
  lines(y = rep(R_stats$specificity,length(seq_)),x =seq_,lty=1,lwd = lwd_,pch=18,col = "darkgray")  
  lines(pears_spec,x =seq_,lty=line_type[2],lwd = lwd_,col = "darkgray")  
  lines(spearman_spec,x =seq_,lty=line_type[3],lwd = lwd_,col = "darkgray")  
  lines(kendall_spec,x =seq_,lty=line_type[4],lwd = lwd_,col = "darkgray")  
  
  lines(y = rep(R_stats$sensitivity,length(seq_)),x =seq_,lty=1,lwd = lwd_)  
  lines(pears_sens,x =seq_,lty=line_type[2],lwd = lwd_)  
  lines(spearman_sens,x =seq_,lty=line_type[3],lwd = lwd_)  
  lines(kendall_sens,x =seq_,lty=line_type[4],lwd = lwd_)  
  legend(x = "bottomright",ncol=1,bty="n",
         legend = c("Sensitivity","Specificity","RSC","Pearson","Spaerman","Kendall"),
         cex = 0.8,lty=c(NA,NA,colors),fill=c("black","darkgray",NA,NA,NA,NA),border = FALSE)
  axis(1, at=seq_, labels=seq_)
  
  
}