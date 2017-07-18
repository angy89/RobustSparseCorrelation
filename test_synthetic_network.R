source("cormad.R")
source("bin_matrix.R")
source("thresholding.R")
source("compute_statistics.R")
source("correlations_cv.R")


library(parallel)
library(foreach)
library(doMC)
library(grndata)

data = syntren300.data
data = as.matrix(data)
gold = syntren300.net
gold = as.matrix(gold)

n_cores = 32
nTest = 1
samples_sizes = 800
nsplits = 100

for(sample_size in samples_sizes){
  dir.create(path = paste("dati/sample_size",sample_size,sep=""))
  print(paste("Sample size: ",sample_size,"\n"))

  M_list = list()
  Th_mat = matrix(0,5,nTest)
  TBW_na = c()
  rownames(Th_mat) = c("RSC","TBW","Pearson","Spearman","Kendall")
  
  for(index in 1:nTest){
    chips = sample(x = 1:nrow(data),size = sample_size,replace = FALSE)
    red_data = data[chips,]
    
    ################## CORMAD CORRELATION THRESHOLDING 
    print("Evaluating RSC Correlation...\n")
    start.time <- Sys.time()
    RSC = thresholding(data = red_data,nsplits = nsplits,n_cores = n_cores,monitor = 0)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(paste("Time tanken for thresholding: ",time.taken))
    Th_mat["RSC",index] = RSC$threshold
    R_vec = cormad.vec(red_data,jobs = n_cores)
    R = cormad.vec2mat(R_vec)
    R_bin = bin_matrix(R,RSC$threshold)
       
    ################## PEARSON CORRELATION THRESHOLDING 
    print("Evaluating Pearson Correlation...\n")
    start.time <- Sys.time()
    pears_th = cv.Rn(data = red_data,nsplits = nsplits,monitor = 0,n_cores = n_cores)
    Th_mat["Pearson",index] = pears_th$threshold
    time.taken <- end.time - start.time
    print(paste("Time tanken for pearson correlation thresholding: ",time.taken))
    print("Evaluating Pearson Correlation...\n")
    pears = cor(red_data,method = "pearson")
    pears_bin =  bin_matrix(pears,pears_th$threshold)
    
    ######################### SPEARMAN CORRELATION THRESHOLDIGN
    print("Evaluating Spearman Correlation...\n")
    start.time <- Sys.time()
    spearman_th = cv.Rho(data = red_data,nsplits = nsplits,monitor = 0,n_cores = n_cores)
    Th_mat["Spearman",index] = spearman_th$threshold
    time.taken <- end.time - start.time
    print(paste("Time tanken for spearman correlation thresholding: ",time.taken))
    print("Evaluating Spearman Correlation...\n")
    spearman = cor(red_data,method = "spearman")
    spearman_bin = bin_matrix(spearman,spearman_th$threshold)

    ######################### KENDALL CORRELATION THRESHOLDIGN
    print("Evaluating Kendall Correlation...\n")
    start.time <- Sys.time()
    kendall_th = cv.Tau(data = red_data,nsplits = nsplits,monitor = 0,n_cores = n_cores)
    Th_mat["Kendall",index] = kendall_th$threshold
    time.taken <- end.time - start.time
    print(paste("Time tanken for kendall correlation thresholding: ",time.taken))
    print("Evaluating Kendall Correlation...\n")
    kendall = cor(red_data,method = "kendall")
    kendall_bin = bin_matrix(kendall,kendall_th$threshold)
    
    ######################### CREATE LIST OF CORRELATION MATRICES
    corList = list(R_bin,pears_bin,spearman_bin,kendall_bin)
    names(corList) = c("RSC","Pearson","Spearman","Kendall")
    
    #########################  COMPARE STATISTICS
    #compare_statistics(R_bin,pears,spearman,kendall,gold)
    M = compare_statistics_after_threshold(listMat = corList,gold)
    

  }

}#end for sample size






