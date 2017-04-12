library("parallel")
source("cormad.R")
source("thresholding.R")
source("compute_statistics.R")
source("find_indices_in_cor_and_true_net.R")

library("grndata")
data = syntren300.data
gold = syntren300.net
colSumsFlosses = c()

sample_size = 160
nsplits = 10
n_cores = 8

chips = sample(x = 1:nrow(data),size = sample_size,replace = FALSE)
red_data = data[chips,]

FLOSSES = thresholding(data = red_data,nsplits = nsplits,n_cores = n_cores)
colSumsFlosses = rbind(colSumsFlosses,colSums(FLOSSES))
th = which.min(colSums(FLOSSES))/100

print("Evaluating Robust Correlation...\n")
R_vec = cormad.vec(red_data,jobs = n_cores)
R = cormad.vec2mat(R_vec)

R_bin = R
R_bin[abs(R)<=th] = 0
R_bin[abs(R)>th] = 1

print("Evaluating Pearson Correlation...\n")
pears = cor(data[chips,],method = "pearson")

print("Evaluating Spearman Correlation...\n")
spearman = cor(data[chips,],method = "spearman")

print("Evaluating Kendall Correlation...\n")
kendall = cor(data,method = "kendall")

corList = list(R,pears,spearman,kendall)

compare_statistics(R_bin,pears,spearman,kendall,gold)

par(mfrow=c(2,2))
robcor = find_indices_in_cor_and_true_net(corList[[1]],gold,th=th,main_text = "Robust")
robcor = find_indices_in_cor_and_true_net(corList[[2]],gold,th=th,main_text = "Pearson")
robcor = find_indices_in_cor_and_true_net(corList[[3]],gold,th=th,main_text = "Spearman")
robcor = find_indices_in_cor_and_true_net(corList[[4]],gold,th=th,main_text = "Kendall")

  