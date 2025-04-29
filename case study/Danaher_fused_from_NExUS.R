rm(list=ls())
set.seed(1)
#install.packages("JGL")
#install.packages("mvtnorm")
#install.packages("psych")
#install.packages("doParallel")
#install.packages("parallel")
library(JGL)
library(mvtnorm)
library(psych)
library(doParallel)
library(parallel)
no_cores <- detectCores()

setwd("Y:/MATLAB codes/Graph Lasso/Jan_12_onwards_size_40_03_07_beta")
no_of_simulations <- 100
C <- 4
p <- 20
Precisions <- list()
for(c in 1:C)
{file <- paste0("A",c,"_sim1.csv")
Precisions[[c]] <- read.csv(file, header = F)
}
sample_sizes <- c(20,40,60,80)
ALL_samples <- list()


for(seed_no in 1:no_of_simulations)
{sample <- list()
for(c in 1:C)
{sample[[c]] <- rmvnorm(sample_sizes[c], mean = rep(0,p), sigma = solve(Precisions[[c]]))}
ALL_samples[[seed_no]] <- sample}



find_opt_lambda <- function(seed_no)
{library(JGL)
  library(psych)
  sample <- ALL_samples[[seed_no]]
  total_lambda_divisions <- 50
  start_lambda <- 0.01
  end_lambda <- 1
  lambda_ones <- exp(seq(from = log(start_lambda), to = log(end_lambda), length.out= total_lambda_divisions))
  lambda_twos <- exp(seq(from = log(start_lambda), to = log(end_lambda), length.out = total_lambda_divisions))
  AIC_all_possible <- matrix(99999999999,total_lambda_divisions,total_lambda_divisions)
  for(i in 1:total_lambda_divisions)
  {for(j in 1:total_lambda_divisions)
  {print(c(seed_no,i,j))
    sprintf('tasks completed: %d; tag: %d\n', i, j)
    JGL_result <- JGL(sample,penalty="fused",lambda1 = lambda_ones[i],lambda2=lambda_twos[j])
  Precision_estimated <- JGL_result$theta
  if(length(Precision_estimated[[1]]) == p^2)
  {num_non_zero_edges <- rep(0,C)
  AIC_sum <- 0
  for(c in 1:C)
  {#print(c(i,j,c))
    num_non_zero_edges[c] <- (length(which(abs(Precision_estimated[[c]])>0.0001))-20)/2
    S_mat <- t(sample[[c]])%*%sample[[c]]
    AIC_sum <- AIC_sum + sample_sizes[c]*tr(S_mat%*%Precision_estimated[[c]]) - 
      sample_sizes[c]*log(det(Precision_estimated[[c]])) + 2*num_non_zero_edges[c]}
  AIC_all_possible[i,j] <- AIC_sum}
  }
  }
  #AIC_all_possible
  min_index <- which(AIC_all_possible == min(AIC_all_possible), arr.ind = TRUE)
  lambda_opt <- c(lambda_ones[min_index[1]],lambda_twos[min_index[2]])
  return(list(lambda_opt = lambda_opt, AICs = AIC_all_possible))}





cl <- makeCluster(no_cores,outfile="")
registerDoParallel(cl)
ptm <- proc.time()
writeLines(c(""), "log.txt")
ANSWER <- foreach(seed_no = 1:no_of_simulations) %dopar% {
  sink("log.txt", append=TRUE)
  find_opt_lambda(seed_no)}
proc.time() - ptm



ALL_lambda <- matrix(0, no_of_simulations,2)
for(i in 1:no_of_simulations)
{ALL_lambda[i,] <- ANSWER[[i]]$lambda_opt}



write.table(ALL_lambda, file = "LAMBDAS_fusedGL.csv", sep=",", row.names = FALSE, col.names=FALSE)
#####################################################################################

ALL_lambda <- read.csv("LAMBDAS_fusedGL.csv",  header = F)


lambda_opt_array <- matrix(0,no_of_simulations,2)
Precision_estimated_array <- list()
for(seed_no in 1:no_of_simulations)
{JGL_result_final <- JGL(ALL_samples[[seed_no]],penalty="fused",lambda1 = ALL_lambda[seed_no,1],
                         lambda2=ALL_lambda[seed_no,2])
Precision_estimated_array[[seed_no]] <- JGL_result_final$theta}



############ True 0-1 array ################
zero_one_array_true <- list()
diag_serial_nos <-  seq(from = 1, to = p^2, by = (p+1))
for(c in 1:C)
{zero_one_array_true[[c]] <- as.vector(t(Precisions[[c]]))
for(i in 1:(p^2))
{if(abs(zero_one_array_true[[c]][i]) > 0)
{zero_one_array_true[[c]][i] <- 1}
  else
  {zero_one_array_true[[c]][i] <- 0}
}
zero_one_array_true[[c]] <- zero_one_array_true[[c]][-diag_serial_nos]
}

intersect_cell_ones <- as.list(numeric(C^2))
intersect_cell_zeros <- as.list(numeric(C^2))
dim(intersect_cell_ones) <- c(C,C)
dim(intersect_cell_zeros) <- c(C,C)


for(c in 1:C)
{for(c_prime in 1:C)
{intersect_cell_ones[[c,c_prime]] <- intersect(which(zero_one_array_true[[c]] == 1),which(zero_one_array_true[[c_prime]] == 1))
intersect_cell_zeros[[c,c_prime]] <- intersect(which(zero_one_array_true[[c]] == 0),which(zero_one_array_true[[c_prime]] == 0))
}}



########### MCC ##################################
zero_one_array <- list()
TPR <- FPR <- TNR <- FNR <- TP <- FP <- TN <- FN <- list()
TPR_COMMON <- FPR_COMMON <- TP_COMMON <- FP_COMMON <- TN_COMMON <- FN_COMMON <- list()
MCC <- list()
MCC_COMMON <- list()


for(seed_no in 1:no_of_simulations)
{for(c in 1:C)
{temp <- as.vector(t(Precision_estimated_array[[seed_no]][[c]]))
zero_one_array[[c]] <- rep(0,p^2)
for(i in 1:(p^2))
{if(abs(temp[i])>0.00001)
{zero_one_array[[c]][i] <- 1}
  else
  {zero_one_array[[c]][i] <- 0}
}
zero_one_array[[c]] <- zero_one_array[[c]][-diag_serial_nos]
}
  
  
  TPR[[seed_no]] <- rep(0,C)
  FPR[[seed_no]] <- rep(0,C)
  TNR[[seed_no]] <- rep(0,C)
  FNR[[seed_no]] <- rep(0,C)
  TP[[seed_no]] <- rep(0,C)
  FP[[seed_no]] <- rep(0,C)
  TN[[seed_no]] <- rep(0,C)
  FN[[seed_no]] <- rep(0,C)
  
  
  true_positive_sum <- rep(0,C)
  true_zero_sum <- rep(0,C)
  estimated_true_positive_sum <- rep(0,C)
  estimated_false_negative_sum <- rep(0,C)
  estimated_false_positive_sum <- rep(0,C)
  estimated_true_negative_sum <- rep(0,C)
  
  for(c in 1:C)
  {for(i in 1:(p^2-p))
  {if(zero_one_array_true[[c]][i] == 1)
  {true_positive_sum[c] <- true_positive_sum[c] + 1
  if(zero_one_array[[c]][i] == 1)
  {estimated_true_positive_sum[c] <- estimated_true_positive_sum[c] + 1}
  if(zero_one_array[[c]][i] == 0)
  {estimated_false_negative_sum[c] <- estimated_false_negative_sum[c] + 1}
  }
  }
    for(i in 1:(p^2-p))
    {if(zero_one_array_true[[c]][i] == 0)
    {true_zero_sum[c] = true_zero_sum[c] + 1
    if(zero_one_array[[c]][i] == 1)
    {estimated_false_positive_sum[c] <- estimated_false_positive_sum[c] + 1}
    if(zero_one_array[[c]][i] == 0)
    {estimated_true_negative_sum[c] <- estimated_true_negative_sum[c] + 1}
    }
    }
    TPR[[seed_no]][c] <- estimated_true_positive_sum[c]/true_positive_sum[c]
    FNR[[seed_no]][c] <- estimated_false_negative_sum[c]/true_positive_sum[c]
    FPR[[seed_no]][c] <- estimated_false_positive_sum[c]/true_zero_sum[c]
    TNR[[seed_no]][c] <- estimated_true_negative_sum[c]/true_zero_sum[c]
    
    TP[[seed_no]][c] <- estimated_true_positive_sum[c]
    FN[[seed_no]][c] <- estimated_false_negative_sum[c]
    FP[[seed_no]][c] <- estimated_false_positive_sum[c]
    TN[[seed_no]][c] <- estimated_true_negative_sum[c]
  }
  
  
  
  TPR_COMMON[[seed_no]] <- matrix(0,C,C)
  FPR_COMMON[[seed_no]] <- matrix(0,C,C)
  
  TP_COMMON[[seed_no]] <- matrix(0,C,C)
  FP_COMMON[[seed_no]] <- matrix(0,C,C)
  TN_COMMON[[seed_no]] <- matrix(0,C,C)
  FN_COMMON[[seed_no]] <- matrix(0,C,C)
  
  for(c in 1:(C-1))
  {for(c_prime in (c+1):C)
  {length_now <- length(intersect_cell_ones[[c,c_prime]])
  if(length_now > 0)
  {TP_COMMON[[seed_no]][c,c_prime] <- sum(zero_one_array[[c]][intersect_cell_ones[[c,c_prime]]]+
                                            zero_one_array[[c_prime]][intersect_cell_ones[[c,c_prime]]])
  FN_COMMON[[seed_no]][c,c_prime] <- 2*length_now - TP_COMMON[[seed_no]][c,c_prime]
  TPR_COMMON[[seed_no]][c,c_prime] <- TP_COMMON[[seed_no]][c,c_prime]/(2*length_now)}
  else
  {TP_COMMON[[seed_no]][c,c_prime] <- NA
  FN_COMMON[[seed_no]][c,c_prime] <- NA
  TPR_COMMON[[seed_no]][c,c_prime] <- NA}
  
  length_now_2 <- length(intersect_cell_zeros[[c,c_prime]])
  if(length_now_2 > 0)
  {FP_COMMON[[seed_no]][c,c_prime] <- sum(zero_one_array[[c]][intersect_cell_zeros[[c,c_prime]]]+
                                            zero_one_array[[c_prime]][intersect_cell_zeros[[c,c_prime]]])
  TN_COMMON[[seed_no]][c,c_prime] <- 2*length_now_2 - FP_COMMON[[seed_no]][c,c_prime]
  FPR_COMMON[[seed_no]][c,c_prime] <- FP_COMMON[[seed_no]][c,c_prime]/(2*length_now_2)}
  else
  {FP_COMMON[[seed_no]][c,c_prime] <- NA
  TN_COMMON[[seed_no]][c,c_prime] <- NA
  FPR_COMMON[[seed_no]][c,c_prime] <- NA}
  }}
  
  
  
  MCC[[seed_no]] <- rep(0,C)
  
  for(i in 1:C)
  {MCC[[seed_no]][i] <- (TP[[seed_no]][i]*TN[[seed_no]][i] - 
                           FP[[seed_no]][i]*FN[[seed_no]][i])/sqrt((TP[[seed_no]][i] + FP[[seed_no]][i])*
                                                                     (TP[[seed_no]][i] + FN[[seed_no]][i])*
                                                                     (TN[[seed_no]][i] + FP[[seed_no]][i])*
                                                                     (TN[[seed_no]][i] + FN[[seed_no]][i]))
  }
  
  
  MCC_COMMON[[seed_no]] <- matrix(0,C,C)
  
  for(c in 1:(C-1))
  {for(c_prime in (c+1):C)
  {MCC_COMMON[[seed_no]][c,c_prime] <- (TP_COMMON[[seed_no]][c,c_prime]*TN_COMMON[[seed_no]][c,c_prime] - 
                                          FP_COMMON[[seed_no]][c,c_prime]*FN_COMMON[[seed_no]][c,c_prime])/
    sqrt((TP_COMMON[[seed_no]][c,c_prime] + FP_COMMON[[seed_no]][c,c_prime])*
           (TP_COMMON[[seed_no]][c,c_prime] + FN_COMMON[[seed_no]][c,c_prime])*
           (TN_COMMON[[seed_no]][c,c_prime] + FP_COMMON[[seed_no]][c,c_prime])*
           (TN_COMMON[[seed_no]][c,c_prime] + FN_COMMON[[seed_no]][c,c_prime]))
  }}
  
}


MCC_mean <- rep(0,C)
MCC_sd <- rep(0,C)

for(i in 1:C)
{temp_array <- rep(0,no_of_simulations)
for(seed_no in 1:no_of_simulations)
{temp_array[seed_no] <- MCC[[seed_no]][i]}
NAs <- which(is.na(temp_array))
if(length(NAs)>0)
{temp_array <- temp_array[-which(is.na(temp_array))]}
MCC_mean[i] <- mean(temp_array)
MCC_sd[i] <- sd(temp_array)
}




MCC_COMMON_mean_structured <- matrix(0,C,C)
MCC_COMMON_sd_structured <- matrix(0,C,C)

for(c in 1:(C-1))
{for(c_prime in (c+1):C)
{temp_array <- rep(0,no_of_simulations)
for(seed_no in 1:no_of_simulations)
{temp_array[seed_no] <- MCC_COMMON[[seed_no]][c,c_prime]}
NAs <- which(is.na(temp_array))
if(length(NAs)>0)
{temp_array <- temp_array[-which(is.na(temp_array))]}
MCC_COMMON_mean_structured[c,c_prime] <- mean(temp_array)
MCC_COMMON_sd_structured[c,c_prime] <- sd(temp_array)
}}

MCC_mean <- round(100*MCC_mean)/100
MCC_sd <- round(100*MCC_sd)/100
MCC_COMMON_mean_structured <- round(100*MCC_COMMON_mean_structured)/100
MCC_COMMON_sd_structured <- round(100*MCC_COMMON_sd_structured)/100


MCC_COMMON_mean_array <- rep(0,C*(C-1)/2)
MCC_COMMON_sd_array <- rep(0,C*(C-1)/2)

count <- 0
for(c in 1:(C-1))
{for(c_prime in (c+1):C)
{count <- count + 1
MCC_COMMON_mean_array[count] <- MCC_COMMON_mean_structured[c,c_prime]
MCC_COMMON_sd_array[count] <- MCC_COMMON_sd_structured[c,c_prime]}}


ALL_MCC_mean <- cbind(t(MCC_mean),t(MCC_COMMON_mean_array))
ALL_MCC_sd <- cbind(t(MCC_sd),t(MCC_COMMON_sd_array))

ALL_MCC_mean_sd <- rbind(ALL_MCC_mean,ALL_MCC_sd)


write.table(ALL_MCC_mean_sd, file = "ZZZZ_MCC_ALL_fusedGL.csv", sep = ",",row.names = FALSE,col.names=FALSE)




# system.time(save1 <- lapply(1:4, find_opt_lambda))
# system.time(save2 <- mclapply(1:4, find_opt_lambda))



