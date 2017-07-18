source("R/AdaptiveThresholding/cormad.R")
source("R/AdaptiveThresholding/bin_matrix.R")
source("R/AdaptiveThresholding/thresholding.R")
source("R/AdaptiveThresholding/compute_statistics.R")
source("R/AdaptiveThresholding/correlations_cv.R")


library(parallel)
library(foreach)
library(doMC)
library(grndata)

data = syntren300.data
data = as.matrix(data)
gold = syntren300.net
gold = as.matrix(gold)

n_cores = 32
samples_sizes = 60
nsplits = 100

chips = sample(x = 1:nrow(data),size = sample_size,replace = FALSE)
red_data = data[chips,]

################## CORMAD CORRELATION THRESHOLDING 
print("Evaluating RSC Correlation...\n")
RSC = thresholding(data = red_data,nsplits = nsplits,n_cores = n_cores,monitor = 0)
R_vec = cormad.vec(red_data,jobs = n_cores)
R = cormad.vec2mat(R_vec)
R_bin = bin_matrix(R,RSC$threshold)

################## PEARSON CORRELATION THRESHOLDING 
print("Evaluating Pearson Correlation...\n")
pears_th = cv.Rn(data = red_data,nsplits = nsplits,monitor = 0,n_cores = n_cores)
pears = cor(red_data,method = "pearson")
pears_bin =  bin_matrix(pears,pears_th$threshold)

######################### SPEARMAN CORRELATION THRESHOLDIGN
print("Evaluating Spearman Correlation...\n")
spearman_th = cv.Rho(data = red_data,nsplits = nsplits,monitor = 0,n_cores = n_cores)
spearman = cor(red_data,method = "spearman")
spearman_bin = bin_matrix(spearman,spearman_th$threshold)

######################### KENDALL CORRELATION THRESHOLDIGN
print("Evaluating Kendall Correlation...\n")
kendall_th = cv.Tau(data = red_data,nsplits = nsplits,monitor = 0,n_cores = n_cores)
kendall = cor(red_data,method = "kendall")
kendall_bin = bin_matrix(kendall,kendall_th$threshold)

######################### CREATE LIST OF CORRELATION MATRICES
corList = list(R_bin,pears_bin,spearman_bin,kendall_bin)
names(corList) = c("RSC","Pearson","Spearman","Kendall")

#########################  COMPARE STATISTICS
#compare_statistics(R_bin,pears,spearman,kendall,gold)
M = compare_statistics_after_threshold(listMat = corList,gold)



