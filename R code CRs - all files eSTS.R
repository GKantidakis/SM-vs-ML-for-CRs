
# Machine learning versus statistical models for competing risks: development and
# validation of prognostic models in extremity soft-tissue sarcoma (eSTS)

# R code


######################################################################
######################################################################
# file "functions_nncr_eSTSmany.R" (necessary to call it when tuning
# the ML techniques (training data) or assessing the performance
# of the models (validation data)

# contains functions for for data pre-processing (for PLANNCR)
# and the estimation of the predictive performance 
# (discrimination, calibration) for all methods. 
######################################################################
######################################################################

## start of the file

library(survcomp)
library(fastDummies)
library(mstate)


# To install packages
# install_packages <- c("fastDummies", "mstate")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("survcomp")


###############################################################

# function that creates data in the long format for training set (PLANNCR original / extended)

train_creator_eSTSmany <- function(data){
  
  N <- nrow(data)
  # assign survival times to 11 intervals, each is a 1-year period (the last 10+ years)
  data$interval <- cut(data$time, breaks = c(0:10, max(data$time, 11)), labels = FALSE)
  data$survival <- cut(data$time, breaks = c(0:10, max(data$time, 11)), labels = FALSE)
  n.times <- data$interval
  data_long <- data[rep(seq_len(N), times = n.times), ]
  
  # create the correct interval enumeration
  for(i in unique(data_long$id)) {
    n_length <- length(data_long$interval[data_long$id == i])
    data_long$interval[data_long$id == i] <- 1:n_length
  }
  
  data_long$status_long <- vector(mode = "numeric",
                                  length = nrow(data_long))
  
  # put indication 1 on status at the interval that patient has relapsed from eSTSmany
  # and 2 that died
  for (i in 1:nrow(data_long)) {
    if (data_long$status[i] != data_long$status_long[i] &&
        data_long$survival[i] == data_long$interval[i])
      data_long$status_long[i] <- data_long$status[i]
  }
  
  events <- dummy_cols(as.factor(data_long$status_long))
  colnames(events) <- gsub(".data", "event", colnames(events))
  data_long <- data.frame(data_long, events[, seq(2, ncol(events), by = 1)])
  
  # interval scaled
  data_long$interval_scaled <- scale(data_long$interval)
  
  data_long <- data_long[, c("id", "time", "status", "sex", "age", "size", "margin",
                             "chemo", "grade", "histology", "depth", "radioth",
                             "interval", "interval_scaled", "survival",
                             "status_long", "event_0", "event_1", "event_2" )]
  return(data_long)
  
}


###############################################################

# function that creates data in the long format for test set (PLANNCR original / extended)
test_creator_eSTSmany <- function(data){
  
  N <- nrow(data)
  # assign survival times to 11 intervals (the last 10+ years)
  data$interval <- max(cut(data$time,
                           breaks = c(0:10, max(data$time, 11)),
                           labels = FALSE)) # replicate the test data 
  
  
  
  # the true interval survival
  data$survival <- cut(data$time,
                       breaks = c(0:10, max(data$time, 11)),
                       labels = FALSE) 
  
  
  data$id <- 1001:(1000 + N) # define the patient ids abstractly
  n.times <- data$interval
  data_long <- data[rep(seq_len(N), times = n.times), ]
  
  # create the correct interval enumeration
  for(i in unique(data_long$id)) {
    n_length <- length(data_long$interval[data_long$id == i])
    data_long$interval[data_long$id == i] <- 1:n_length
    
  }
  
  data_long$status_long <- vector(mode = "numeric",
                                  length = nrow(data_long))
  
  # put indication 1 on status at the intervals on
  # which a patient has relapsed from eSTSmany and 2 that died
  for (i in 1:nrow(data_long)) {
    if (data_long$status[i] != data_long$status_long[i] && 
        data_long$survival[i] <= data_long$interval[i]) 
      data_long$status_long[i] <- data_long$status[i]
  }
  
  
  events <- dummy_cols(as.factor(data_long$status_long))
  colnames(events) <- gsub(".data", "event", colnames(events))
  data_long <- data.frame(data_long, events[, seq(2, ncol(events), by = 1)])
  
  # interval scaled
  data_long$interval_scaled <- scale(data_long$interval)
  
  data_long <- data_long[, c("id", "time", "status", "sex", "age", "size", "margin",
                             "chemo", "grade", "histology", "depth", "radioth",
                             "interval", "interval_scaled", "survival",
                             "status_long", "event_0", "event_1", "event_2")]
  return(data_long)
  
}


# function that calculates Brier score and AUC for cause 1

metrics_cause1 <- function(risk_matrix, data) {
  
  f_nn <- Hist(time, status) ~ 1
  # we subtract 1 day at 11 years to avoid zero censoring prob weight at this time  
  bsc <- riskRegression:::Score(object = list("nnet" = risk_matrix),
                                formula = f_nn,
                                data = data,
                                metrics = c("Brier"),
                                cause = 1,
                                times = c(0:10, 10.999),
                                null.model = FALSE,
                                conservative = TRUE, # speed up computation
                                conf.int = FALSE,
                                split.method = "none")
  
  auc <- riskRegression:::Score(object = list("nnet" = risk_matrix),
                                formula = f_nn,
                                data = data,
                                metrics = c("AUC"),
                                cause = 1, 
                                times = c(0:10, 10.999),
                                null.model = FALSE,
                                conservative = TRUE, # speed up computation
                                conf.int = FALSE,
                                split.method = "none")
  auc$AUC$score[, 3][1] <- 1 # assign 1 if NaNs
  
  
  return(list(Time = seq(0, 11, length.out = 12), 
              Brier = as.numeric(unlist(bsc$Brier$score[, 3])),
              AUC = as.numeric(unlist(auc$AUC$score[, 3]))))
}


# function that calculates Brier score and AUC for cause 2

metrics_cause2 <- function(risk_matrix, data) {
  
  f_nn <- Hist(time, status) ~ 1
  # we subtract 1 day at 11 years to avoid zero censoring prob weight at this time  
  bsc <- riskRegression:::Score(object = list("nnet" = risk_matrix),
                                formula = f_nn,
                                data = data,
                                metrics = c("Brier"),
                                cause = 2,
                                times = c(0:10, 10.999),
                                null.model = FALSE,
                                conservative = TRUE, # speed up computation
                                conf.int = FALSE,
                                split.method = "none")
  
  auc <- riskRegression:::Score(object = list("nnet" = risk_matrix),
                                formula = f_nn,
                                data = data,
                                metrics = c("AUC"),
                                cause = 2, 
                                times = c(0:10, 10.999),
                                null.model = FALSE,
                                conservative = TRUE, # speed up computation
                                conf.int = FALSE,
                                split.method = "none")
  auc$AUC$score[, 3][1] <- 1 # assign 1 if NaNs
  
  
  return(list(Time = seq(0, 11, length.out = 12), 
              Brier = as.numeric(unlist(bsc$Brier$score[, 3])),
              AUC = as.numeric(unlist(auc$AUC$score[, 3]))))
}


# miscalibration as the MSE of 4 groups (observed minus the predicted cumulative incidence
# for a cause)

calibration_mse <- function(cumprob_vec, data, time, cutpoints = 4, cause) { 
  
  
  if (var(cumprob_vec, na.rm = TRUE) != 0) { # run if not all probabilities are equal
    
    cuts <- unique(quantile(cumprob_vec, seq(0, 1, length = cutpoints + 1), na.rm = TRUE))
    
    if ((length(cuts) - 1) == cutpoints) { # if cutpoints groups do exist
      
      data$groups <- as.factor(as.numeric(cut(x = cumprob_vec, breaks = as.vector(cuts),
                                              include.lowest = TRUE)))
      
      # average cumulative incidence of event c
      average_cumprobs <- vector(mode = "numeric", length = length(cuts) - 1)
      for (i in 1:(length(cuts) - 1)) {
        average_cumprobs[i] <- mean(cumprob_vec[which(as.numeric(data$groups) == i)], na.rm = TRUE)
      }
      
      # compare with observed cumulative incidence of event of type cause
      observed_cumprobs_all <- summary(Cuminc("time", "status", group = "groups", data = data),
                                       times = time, extend = TRUE)
      
      observed_cumprobs <-  observed_cumprobs_all$pstate[, 1 + cause]
      
      
      mse <- sum((observed_cumprobs - average_cumprobs)^2, na.rm = TRUE) / length(average_cumprobs)
      
    } else { # else if they are equal set to NA
      
      mse <- NA
      
    } 
    
  } else { # else if they are equal set to NA
    
    mse <- NA
    
  } 
  
  # return the mean squared errors calculated
  return(list(mse_groups = mse))
  
}



# calculate the predictive performance measures for PLANNCR original

measures_nnet_eSTSmany <- function(trained_model, datashort, datalong) {
  
  test_x <- datalong[, c(4:17, 19)]
  dimnames(test_x) <- NULL
  datalong <- as.data.frame(datalong)
  real_times <- datashort$time
  real_status <- datashort$status
  
  df1 <- data.frame(predict(trained_model, test_x, type = "raw"))
  colnames(df1) <- c("prob_0", "prob_1", "prob_2")
  
  df1$id <- as.vector(datalong$id) # ids of the patients
  df1$interval <- as.vector(datalong$interval)
  df1$survival <- datalong$survival
  groups <- split(df1, f = df1$id) # split to lists per id
  
  # discrete cause-specific cumulative incidence: cause 1
  risks1_probs <-  lapply(groups, function(x) {
    x <- cumsum(cumprod(1 - (x$prob_1 + x$prob_2))*x$prob_1)})
  # cumulative incidence of an event of cause 1 until time t 
  # in the rows the values should be increasing
  risks1_mat <- cbind(0, do.call("rbind", risks1_probs))
  dimnames(risks1_mat) <- NULL
  
  # discrete cause-specific cumulative incidence: cause 2
  risks2_probs <-  lapply(groups, function(x) {
    x <- cumsum(cumprod(1 - (x$prob_1 + x$prob_2))*x$prob_2)})
  risks2_mat <- cbind(0, do.call("rbind", risks2_probs))
  dimnames(risks2_mat) <- NULL
  
  df2 <- data.frame(time = real_times,
                    status = real_status)
  
  # for cause-specific hazards 1
  metrics_obj1 <- metrics_cause1(risk_matrix = risks1_mat, data = df2)
  brier_set1 <- metrics_obj1$Brier
  auc_set1 <- metrics_obj1$AUC
  mse_set1_time2 <- calibration_mse(cumprob_vec = risks1_mat[, 3], data = df2,
                                    time = 2, cutpoints = 4, cause = 1)$mse_groups
  mse_set1_time5 <- calibration_mse(cumprob_vec = risks1_mat[, 6], data = df2,
                                    time = 5, cutpoints = 4, cause = 1)$mse_groups
  mse_set1_time10 <- calibration_mse(cumprob_vec = risks1_mat[, 11], data = df2,
                                     time = 10, cutpoints = 4, cause = 1)$mse_groups
  
  
  
  # for cause-specific hazards 2
  metrics_obj2 <- metrics_cause2(risk_matrix = risks2_mat, data = df2)
  brier_set2 <- metrics_obj2$Brier
  auc_set2 <- metrics_obj2$AUC
  mse_set2_time2 <- calibration_mse(cumprob_vec = risks2_mat[, 3], data = df2,
                                    time = 2, cutpoints = 4, cause = 2)$mse_groups
  mse_set2_time5 <- calibration_mse(cumprob_vec = risks2_mat[, 6], data = df2,
                                    time = 5, cutpoints = 4, cause = 2)$mse_groups
  mse_set2_time10 <- calibration_mse(cumprob_vec = risks2_mat[, 11], data = df2,
                                     time = 10, cutpoints = 4, cause = 2)$mse_groups
  
  return(list(Risks1_mat = risks1_mat,
              Risks2_mat = risks2_mat,
              Time = metrics_obj1$Time,
              Brier_cause1 = brier_set1,
              AUC_cause1 = auc_set1,
              MSE2_cause1 = mse_set1_time2,
              MSE5_cause1 = mse_set1_time5,
              MSE10_cause1 = mse_set1_time10,
              Brier_cause2 = brier_set2,
              AUC_cause2 = auc_set2,
              MSE2_cause2 = mse_set2_time2,
              MSE5_cause2 = mse_set2_time5,
              MSE10_cause2 = mse_set2_time10))
  
}


# calculate the predictive performance measures for PLANNCR extended

measures_keras_eSTSmany <- function(trained_model, datashort, datalong) {
  
  test_x <- datalong[, c(4:17, 19)]
  dimnames(test_x) <- NULL
  datalong <- as.data.frame(datalong)
  real_times <- datashort$time
  real_status <- datashort$status
  
  df1 <- data.frame(predict(trained_model, test_x, batch_size = 32))
  colnames(df1) <- c("prob_0", "prob_1", "prob_2")
  
  df1$id <- as.vector(datalong$id) # ids of the patients
  df1$interval <- as.vector(datalong$interval)
  df1$survival <- datalong$survival
  groups <- split(df1, f = df1$id) # split to lists per id
  
  # discrete cause-specific cumulative incidence: cause 1
  risks1_probs <-  lapply(groups, function(x) {
    x <- cumsum(cumprod(1 - (x$prob_1 + x$prob_2))*x$prob_1)})
  # cumulative incidence of an event of cause 1 until time t 
  # in the rows the values should be increasing
  risks1_mat <- cbind(0, do.call("rbind", risks1_probs))
  dimnames(risks1_mat) <- NULL
  
  # discrete cause-specific cumulative incidence: cause 2
  risks2_probs <-  lapply(groups, function(x) {
    x <- cumsum(cumprod(1 - (x$prob_1 + x$prob_2))*x$prob_2)})
  risks2_mat <- cbind(0, do.call("rbind", risks2_probs))
  dimnames(risks2_mat) <- NULL
  
  df2 <- data.frame(time = real_times,
                    status = real_status)
  
  # for cause-specific hazards 1
  metrics_obj1 <- metrics_cause1(risk_matrix = risks1_mat, data = df2)
  brier_set1 <- metrics_obj1$Brier
  auc_set1 <- metrics_obj1$AUC
  mse_set1_time2 <- calibration_mse(cumprob_vec = risks1_mat[, 3], data = df2,
                                    time = 2, cutpoints = 4, cause = 1)$mse_groups
  mse_set1_time5 <- calibration_mse(cumprob_vec = risks1_mat[, 6], data = df2,
                                    time = 5, cutpoints = 4, cause = 1)$mse_groups
  mse_set1_time10 <- calibration_mse(cumprob_vec = risks1_mat[, 11], data = df2,
                                     time = 10, cutpoints = 4, cause = 1)$mse_groups
  
  
  # for cause-specific hazards 2
  metrics_obj2 <- metrics_cause2(risk_matrix = risks2_mat, data = df2)
  brier_set2 <- metrics_obj2$Brier
  auc_set2 <- metrics_obj2$AUC
  mse_set2_time2 <- calibration_mse(cumprob_vec = risks2_mat[, 3], data = df2,
                                    time = 2, cutpoints = 4, cause = 2)$mse_groups
  mse_set2_time5 <- calibration_mse(cumprob_vec = risks2_mat[, 6], data = df2,
                                    time = 5, cutpoints = 4, cause = 2)$mse_groups
  mse_set2_time10 <- calibration_mse(cumprob_vec = risks2_mat[, 11], data = df2,
                                     time = 10, cutpoints = 4, cause = 2)$mse_groups
  
  return(list(Risks1_mat = risks1_mat,
              Risks2_mat = risks2_mat,
              Time = metrics_obj1$Time,
              Brier_cause1 = brier_set1,
              AUC_cause1 = auc_set1,
              MSE2_cause1 = mse_set1_time2,
              MSE5_cause1 = mse_set1_time5,
              MSE10_cause1 = mse_set1_time10,
              Brier_cause2 = brier_set2,
              AUC_cause2 = auc_set2,
              MSE2_cause2 = mse_set2_time2,
              MSE5_cause2 = mse_set2_time5,
              MSE10_cause2 = mse_set2_time10))
  
}


# calculate the predictive performance for cause-specific Cox, Fine-Gray, RSFCR

measures_eSTSmany <- function(risk_mat1, risk_mat2, datashort) {
  
  # for cause-specific hazards 1
  metrics_obj1 <- metrics_cause1(risk_matrix = risk_mat1, data = datashort)
  brier_set1 <- metrics_obj1$Brier
  auc_set1 <- metrics_obj1$AUC
  mse_set1_time2 <- calibration_mse(cumprob_vec = risk_mat1[, 3], data = datashort,
                                    time = 2, cutpoints = 4, cause = 1)$mse_groups
  mse_set1_time5 <- calibration_mse(cumprob_vec = risk_mat1[, 6], data = datashort,
                                    time = 5, cutpoints = 4, cause = 1)$mse_groups
  mse_set1_time10 <- calibration_mse(cumprob_vec = risk_mat1[, 11], data = datashort,
                                     time = 10, cutpoints = 4, cause = 1)$mse_groups
  
  
  # for cause-specific hazards 2
  metrics_obj2 <- metrics_cause2(risk_matrix = risk_mat2, data = datashort)
  brier_set2 <- metrics_obj2$Brier
  auc_set2 <- metrics_obj2$AUC
  mse_set2_time2 <- calibration_mse(cumprob_vec = risk_mat2[, 3], data = datashort,
                                    time = 2, cutpoints = 4, cause = 2)$mse_groups
  mse_set2_time5 <- calibration_mse(cumprob_vec = risk_mat2[, 6], data = datashort,
                                    time = 5, cutpoints = 4, cause = 2)$mse_groups
  mse_set2_time10 <- calibration_mse(cumprob_vec = risk_mat2[, 11], data = datashort,
                                     time = 10, cutpoints = 4, cause = 2)$mse_groups
  
  return(list(Risks1_mat = risk_mat1,
              Risks2_mat = risk_mat2,
              Time = metrics_obj1$Time,
              Brier_cause1 = brier_set1,
              AUC_cause1 = auc_set1,
              MSE2_cause1 = mse_set1_time2,
              MSE5_cause1 = mse_set1_time5,
              MSE10_cause1 = mse_set1_time10,
              Brier_cause2 = brier_set2,
              AUC_cause2 = auc_set2,
              MSE2_cause2 = mse_set2_time2,
              MSE5_cause2 = mse_set2_time5,
              MSE10_cause2 = mse_set2_time10))
  
}  


## end of the file



######################################################################
######################################################################
# file "PLANNCR_tune_eSTSmany.R" 

# Used to tune the parameters for PLANNCR original,
# PLANNCR extended or RSFCR 
######################################################################
######################################################################

## start of the file


# To install packages
# install_packages <- c("riskRegression", "prodlim", "ggplot2", "survival",
#                        "data.table", "nnet", "keras", "tictoc")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }


library(riskRegression)
library(prodlim)
library(ggplot2)
library(survival)
library(data.table)
library(nnet)
library(keras)
library(tictoc)
source("functions_nncr_eSTSmany.R")
load("eSTSmany.Rdata") # call the collected data

# how to set up keras in R
# first install anaconda 3
# in anaconda command prompt use the following:
# pip install tensorflow
# pip install keras

# test by running: import tensorflow as tf

# this is done
# install.packages("reticulate")
# install.packages("tensorflow")
# install.packages("keras")

# library(tensorflow)
# install_tensorflow()
# library(keras)
# install_keras()


# tuning eSTSmany
# with bootstrap, search for hyperparameters afresh for each loop
# supervised learning steps inside the loop
# B = 100
# original data  ---> N = 3826 develop, N = OOB validate
# p = 9 in original format, 14 with model.matrix (expanding factors to dummy variables)


#############################
# 1. For PLANNCR original 
# tuned with Brier score at 5 years, and AUC at 5 years
#############################


set.seed(12345)
n <- nrow(eSTSmany)
B <- 100

# getting the bootstrap samples without loops
# each column is a bootstrap re-sample
boot_samples <- matrix(sample(eSTSmany$id, size = n*B, replace = TRUE),
                       nrow = n, ncol = B)

# create objects to store the statistics
# create a loop to calculate the statistics

node_size <- seq(2, 14, by = 1) # grid of node sizes
weight_decay <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
combis <- expand.grid(node_size, weight_decay) # 91 combinations

# initialize objects
mat_brier <-  matrix(0, nrow = B, ncol = nrow(combis))
mat_auc <- matrix(0, nrow = B, ncol = nrow(combis))
mat_calib <- matrix(0, nrow = B, ncol = nrow(combis))


# run the bootstrap training
tic("Running the bootstrap training")


for (i in 1:B) {
  
  cat("Started iteration i = ", i, "\n")
  # use the proper bootstrap column
  indices <- boot_samples[, i]
  
  eSTSmany2 <- eSTSmany[indices, ]
  eSTSmany2$id <- 1:nrow(eSTSmany2) 
  
  # sample randomly 3/4 of data as train and 1/4 as test data per bootstrap sample
  set.seed(12345)
  vec <- sample(x = eSTSmany2$id, size = round((3/4)*nrow(eSTSmany2))) 
  eSTSmany2_train <- eSTSmany2[eSTSmany2$id %in% vec, ]
  eSTSmany2_test <- eSTSmany2[!(eSTSmany2$id %in% vec), ]
  
  # CRs and time > 10 years else resample rows 
  e <- 0
  while (sum(eSTSmany2_test$status == 2) == 0 | max(eSTSmany2_test$time) < 10)  {
    
    e <- e + 1
    set.seed(e)
    vec <- sample(x = eSTSmany2$id, size = round((3/4)*nrow(eSTSmany2))) 
    eSTSmany2_train <- eSTSmany2[eSTSmany2$id %in% vec, ]
    eSTSmany2_test <- eSTSmany2[!(eSTSmany2$id %in% vec), ]
  }
  
  # scale data  
  eSTSmany2_train$age <- scale(eSTSmany2_train$age)
  eSTSmany2_train$size <- scale(eSTSmany2_train$size)
  eSTSmany2_test$age <- scale(eSTSmany2_test$age)
  eSTSmany2_test$size <- scale(eSTSmany2_test$size)
  
  train_long_eSTSmany <- train_creator_eSTSmany(data = eSTSmany2_train)
  # categorical variables to numerical (reference level redundant)
  train_long_eSTSmany <- model.matrix(~., train_long_eSTSmany)[, -1] 
  
  train_x_eSTSmany <- train_long_eSTSmany[, c(4:17, 19)] # p = 14 variables + time interval
  dimnames(train_x_eSTSmany) <- NULL # the object must have empty dimnames
  train_y_eSTSmany <- train_long_eSTSmany[, (ncol(train_long_eSTSmany) - 2):ncol(train_long_eSTSmany)]
  dimnames(train_y_eSTSmany) <- NULL
  
  test_long_eSTSmany <- test_creator_eSTSmany(data = eSTSmany2_test)
  # categorical variables to numerical (reference level redundant)
  test_long_eSTSmany <- model.matrix(~., test_long_eSTSmany)[, -1]
  
  test_x_eSTSmany <- test_long_eSTSmany[, c(4:17, 19)] # p = 14 variables + time interval
  dimnames(test_x_eSTSmany) <- NULL # the object must have empty dimnames
  test_y_eSTSmany <- test_long_eSTSmany[, (ncol(train_long_eSTSmany) - 2):ncol(train_long_eSTSmany)]
  dimnames(test_y_eSTSmany) <- NULL
  
  for (j in 1:nrow(combis)) { 
    
    cat("Testing combination number:", j, "of repeat", i, "out of", B, "\n")
    cat("calculating for node size:", combis[j, 1], "and weight decay", combis[j, 2], "...", "\n")
    
    # start building the model in nnet
    
    set.seed(12345)
    fit_nnet <- nnet(x = train_x_eSTSmany, y = train_y_eSTSmany,
                     size = combis[j, 1], decay = combis[j, 2], 
                     maxit = 500, trace = FALSE, softmax = TRUE) # runs automatically with entropy = TRUE
    
    
    results1 <- measures_nnet_eSTSmany(trained_model = fit_nnet,
                                       datashort = eSTSmany2_test, 
                                       datalong = test_long_eSTSmany)
    
    mat_brier[i, j] <- results1$Brier_cause1[6]
    mat_auc[i, j] <- results1$AUC_cause1[6]
    mat_calib[i, j] <- results1$MSE5_cause1
    
    cat("The Brier Score at 5 years is:", round(results1$Brier_cause1[6], 3),
        "and the AUC is:", round(results1$AUC_cause1[6], 3), "\n")
    
  }
  
  
}

time_list <- toc()

cat((time_list$toc - time_list$tic) / 3600, "hours elapsed", "\n") 
# 22.83089 hours elapsed 

df_logistic_ovr <- as.data.frame(cbind(node_size = combis[, 1],
                                       weight_decay = combis[, 2],
                                       brier5 = colMeans(mat_brier),
                                       auc5 = colMeans(mat_auc),
                                       calib5 = colMeans(mat_calib,
                                                         na.rm = TRUE))) 

save(df_logistic_ovr, file = "results_nnet_eSTSmany.RData")

ind <- head(order(df_logistic_ovr$brier5, decreasing = FALSE), 3)
df_logistic_ovr[ind, ]

ind2 <- head(order(df_logistic_ovr$auc5, decreasing = TRUE), 3)
df_logistic_ovr[ind2, ]


#################################################
# 2. For PLANNCR extended (example for sigmoid)
# similarly for relu, tanh activation functions
# tuned with Brier score at 5 years, and AUC at 5 years
#################################################



set.seed(12345)
n <- nrow(eSTSmany)
B <- 100

# getting the bootstrap samples without loops
# each column is a bootstrap re-sample
boot_samples <- matrix(sample(eSTSmany$id, size = n*B, replace = TRUE),
                       nrow = n, ncol = B)

# create objects to store the statistics
# create a loop to calculate the statistics

# 252 combis
node_size <- seq(2, 14, by = 2) # grid of node sizes
dropout_rate <- c(0.1, 0.2, 0.4)
lr <- c(0.1, 0.2, 0.4)
momentum <- c(0.8, 0.9)
class_weights <- c(1, 1.25)
combis <- expand.grid(node_size, dropout_rate, lr, momentum, class_weights)

# initialize objects
mat_brier <-  matrix(0, nrow = B, ncol = nrow(combis))
mat_auc <- matrix(0, nrow = B, ncol = nrow(combis))
mat_calib <- matrix(0, nrow = B, ncol = nrow(combis))

# run the bootstrap training
tic("Running the bootstrap training")


for (i in 94:B) {
  
  cat("Started iteration i = ", i, "\n")
  # use the proper bootstrap column
  indices <- boot_samples[, i]
  
  eSTSmany2 <- eSTSmany[indices, ]
  eSTSmany2$id <- 1:nrow(eSTSmany2) 
  
  # sample randomly 3/4 of data as train and 1/4 as test data per bootstrap sample
  set.seed(12345)
  vec <- sample(x = eSTSmany2$id, size = round((3/4)*nrow(eSTSmany2))) 
  eSTSmany2_train <- eSTSmany2[eSTSmany2$id %in% vec, ]
  eSTSmany2_test <- eSTSmany2[!(eSTSmany2$id %in% vec), ]
  
  # CRs  and time > 10 years else resample rows
  e <- 0
  while (sum(eSTSmany2_test$status == 2) == 0 | max(eSTSmany2_test$time) < 10)  {
    
    e <- e + 1
    set.seed(e)
    vec <- sample(x = eSTSmany2$id, size = round((3/4)*nrow(eSTSmany2))) 
    eSTSmany2_train <- eSTSmany2[eSTSmany2$id %in% vec, ]
    eSTSmany2_test <- eSTSmany2[!(eSTSmany2$id %in% vec), ]
  }
  
  # scale data  
  eSTSmany2_train$age <- scale(eSTSmany2_train$age)
  eSTSmany2_train$size <- scale(eSTSmany2_train$size)
  eSTSmany2_test$age <- scale(eSTSmany2_test$age)
  eSTSmany2_test$size <- scale(eSTSmany2_test$size)
  
  train_long_eSTSmany <- train_creator_eSTSmany(data = eSTSmany2_train)
  # categorical variables to numerical (reference level redundant)
  train_long_eSTSmany <- model.matrix(~., train_long_eSTSmany)[, -1] 
  
  train_x_eSTSmany <- train_long_eSTSmany[, c(4:17, 19)] # p = 14 variables + time interval
  dimnames(train_x_eSTSmany) <- NULL # the object must have empty dimnames
  train_y_eSTSmany <- train_long_eSTSmany[, (ncol(train_long_eSTSmany) - 2):ncol(train_long_eSTSmany)]
  dimnames(train_y_eSTSmany) <- NULL
  
  test_long_eSTSmany <- test_creator_eSTSmany(data = eSTSmany2_test)
  # categorical variables to numerical (reference level redundant)
  test_long_eSTSmany <- model.matrix(~., test_long_eSTSmany)[, -1]
  
  test_x_eSTSmany <- test_long_eSTSmany[, c(4:17, 19)] # p = 14 variables + time interval
  dimnames(test_x_eSTSmany) <- NULL # the object must have empty dimnames
  test_y_eSTSmany <- test_long_eSTSmany[, (ncol(train_long_eSTSmany) - 2):ncol(train_long_eSTSmany)]
  dimnames(test_y_eSTSmany) <- NULL
  
  
  for (j in 1:nrow(combis)) { 
    
    cat("Testing combination number:", j, "of repeat", i, "out of", B, "\n")
    cat("calculating for node size:", combis[j, 1], ", dropout rate:", combis[j, 2], "\n",
        "and learning rate", combis[j, 3], "and momentum", combis[j, 4],
        "and weak class weight", combis[j, 5], "...", "\n")
    
    
    # start building the model in keras
    
    # set global default to never show metrics
    options(keras.view_metrics = FALSE)
    k_clear_session() # to avoid clutter from old models / layers in cross validation
    tensorflow::tf$random$set_seed(12345) # for update to TensorFlow 2.0
    # model with one hidden layer
    
    # Initialize a sequential model
    fit_keras1 <- keras_model_sequential() 
    
    # Add layers to the model
    fit_keras1 %>% 
      layer_dense(units = combis[j, 1], activation = 'sigmoid',
                  input_shape = ncol(train_x_eSTSmany)) %>% 
      layer_dropout(rate = combis[j, 2]) %>%
      layer_dense(units = 3, activation = 'softmax')
    
    # Compile the model
    fit_keras1 %>% compile(
      loss = 'categorical_crossentropy', # for multi-class classification problem
      optimizer = optimizer_sgd(learning_rate = combis[j, 3], momentum = combis[j, 4])
      #, metrics = c("accuracy")
    )
    early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 3) 
    
    result_keras <- fit_keras1 %>% fit(
      train_x_eSTSmany,
      train_y_eSTSmany,
      epochs = 20, 
      batch_size = 32,
      validation_data = list(test_x_eSTSmany, test_y_eSTSmany),
      class_weight = list("0" = 1, "1" = combis[j, 5], "2" = combis[j, 5]),
      callbacks = c(early_stopping),
      verbose = 0 # don't display progress bar
    )
    
    
    results2 <- measures_keras_eSTSmany(trained_model = fit_keras1,
                                        datashort = eSTSmany2_test, 
                                        datalong = test_long_eSTSmany)
    
    mat_brier[i, j] <- results2$Brier_cause1[6]
    mat_auc[i, j] <- results2$AUC_cause1[6]
    mat_calib[i, j] <- results2$MSE5_cause1
    
    cat("The Brier Score at 5 years is:", round(results2$Brier_cause1[6], 3),
        "and the AUC is:", round(results2$AUC_cause1[6], 3), "\n")
    
  }
  
  
}

time_list <- toc()

cat((time_list$toc - time_list$tic) / 3600, "hours elapsed", "\n") 
# 42.12284 hours elapsed

df_sigmoid_ovr <- as.data.frame(cbind(node_size = combis[, 1],
                                      dropout_rate = combis[, 2],
                                      lr = combis[, 3],
                                      momentum = combis[, 4],
                                      weak_class = combis[, 5],
                                      brier5 = colMeans(mat_brier),
                                      auc5 = colMeans(mat_auc),
                                      calib5 = colMeans(mat_calib,
                                                        na.rm = TRUE))) 


save(df_sigmoid_ovr, file = "results_keras_eSTSmany.RData")

# best brier at 5 years
ind <- head(order(df_sigmoid_ovr$brier5, decreasing = FALSE), 5)
df_sigmoid_ovr[ind, ]

# best auc at 5 years
ind2 <- head(order(df_sigmoid_ovr$auc5, decreasing = TRUE), 5)
df_sigmoid_ovr[ind2, ]


#################################################
# 3. For RSFCR, logrank split rule
# tuned with forest error E = (1 - C)
#################################################

set.seed(12345)
n <- nrow(eSTSmany)
B <- 100

# getting the bootstrap samples without loops
# each column is a bootstrap re-sample
boot_samples <- matrix(sample(eSTSmany$id, size = n*B, replace = TRUE),
                       nrow = n, ncol = B)

# create objects to store the statistics
# create a loop to calculate the statistics


mtry <- seq(1, 5, by = 1)
nsplit <- c(2) # seq(2, 3, by = 1)
nodesize <- seq(10, 30, by = 4)
# 30 combinations
combis <- expand.grid(mtry, nsplit, nodesize)

# initialize objects
mat_error <-  matrix(0, nrow = B, ncol = nrow(combis))
mat_calib <- matrix(0, nrow = B, ncol = nrow(combis))

# run the bootstrap training
tic("Running the bootstrap training")


for (i in 1:B) {
  
  cat("Started iteration i = ", i, "\n")
  # use the proper bootstrap column
  indices <- boot_samples[, i]
  
  eSTSmany2 <- eSTSmany[indices, ]
  eSTSmany2$id <- 1:nrow(eSTSmany2) 
  
  # sample randomly 3/4 of data as train and 1/4 as test data per bootstrap sample
  set.seed(12345)
  vec <- sample(x = eSTSmany2$id, size = round((3/4)*nrow(eSTSmany2))) 
  eSTSmany2_train <- eSTSmany2[eSTSmany2$id %in% vec, ]
  eSTSmany2_test <- eSTSmany2[!(eSTSmany2$id %in% vec), ]
  
  # CRs and time > 10 years else resample rows 
  e <- 0
  while (sum(eSTSmany2_test$status == 2) == 0 | max(eSTSmany2_test$time) < 10)  {
    
    e <- e + 1
    set.seed(e)
    vec <- sample(x = eSTSmany2$id, size = round((3/4)*nrow(eSTSmany2))) 
    eSTSmany2_train <- eSTSmany2[eSTSmany2$id %in% vec, ]
    eSTSmany2_test <- eSTSmany2[!(eSTSmany2$id %in% vec), ]
  }
  
  
  for (j in 1:nrow(combis)) { 
    
    cat("Testing combination number:", j, "of repeat", i, "out of", B, "\n")
    cat("calculating for mtry:", combis[j, 1], "nsplit", combis[j, 2],
        "and nodesize", combis[j, 3], "...", "\n")
    
    # start building the rsf
    
    fit_rsf <- rfsrc(Surv(time, status) ~ sex + age + size + margin + chemo 
                     + grade + histology + depth + radioth,
                     data = eSTSmany2_train, ntree = 1000, mtry = combis[j, 1],
                     nsplit = combis[j, 2], nodesize = combis[j, 3],
                     splitrule = "logrank", cause = 1, seed = 12345,
                     forest = TRUE, importance = FALSE)
    
    df <- data.frame(time = eSTSmany2_test$time,
                     status = eSTSmany2_test$status)
    
    r1c <- predictRisk(fit_rsf, newdata=as.data.table(eSTSmany2_test), times=0:11, cause=1)
    r2c <- predictRisk(fit_rsf, newdata=as.data.table(eSTSmany2_test), times=0:11, cause=2)
    
    pred_test <- predict(fit_rsf, newdata = eSTSmany2_test)
    results1 <- measures_eSTSmany(risk_mat1 = r1c, risk_mat2 = r2c, datashort = df)
    
    mat_error[i, j] <- pred_test$err.rate[pred_test$ntree] # OOB error for cause 1 (1 - C)
    mat_calib[i, j] <- results1$MSE5_cause1
    
    cat("The prediction error is:", round(pred_test$err.rate[pred_test$ntree], 3), "\n")
    
  }
  
  
}

time_list <- toc()

cat((time_list$toc - time_list$tic) / (60*60), "hours elapsed", "\n") 
# 15.399 hours elapsed 

df_logrank_ovr <- as.data.frame(cbind(mtry = combis[, 1],
                                      nsplit = combis[, 2],
                                      nodesize = combis[, 3],
                                      error = colMeans(mat_error),
                                      calib5 = colMeans(mat_calib))) 

save(df_logrank_ovr, file = "results_rsf_logrank_eSTSmany.RData")

ind <- head(order(df_logrank_ovr$error, decreasing = FALSE), 3)
df_logrank_ovr[ind, ]


## end of the file



######################################################################
######################################################################
# file "valid_eSTSmany.R" 

# Used to validate the predictive performance for all models
######################################################################
######################################################################


## start of the file


# To install packages
# install_packages <- c("riskRegression", "prodlim", "ggplot2", "survival",
#                        "data.table", "nnet", "keras", "tictoc")
# for (i in 1:length (install_packages)){
#   if (!install_packages[i] %in% installed.packages()){
#     install.packages(install_packages[i])
#   }
# }


library(riskRegression)
library(prodlim)
library(ggplot2)
library(survival)
library(data.table)
library(nnet)
library(keras)
library(tictoc)
source("functions_nncr_eSTSmany.R")
load("eSTSmany.Rdata") # call the collected data

# how to set up keras in R
# first install anaconda 3
# in anaconda command prompt use the following:
# pip install tensorflow
# pip install keras

# test by running: import tensorflow as tf

# this is done
# install.packages("reticulate")
# install.packages("tensorflow")
# install.packages("keras")

# library(tensorflow)
# install_tensorflow()
# library(keras)
# install_keras()


set.seed(12345)
n <- nrow(eSTSmany)
B <- 100

# getting the bootstrap samples without loops
# each column is a bootstrap re-sample
boot_samples <- matrix(sample(eSTSmany$id, size = n*B, replace = TRUE),
                       nrow = n, ncol = B)


# create objects to store the statistics
# create a loop to calculate the statistics

# initialize objects
list_nnet <- vector(mode = "list", length = ncol(boot_samples))  
list_keras <- vector(mode = "list", length = ncol(boot_samples))  
list_csc <- vector(mode = "list", length = ncol(boot_samples))  
list_fg <- vector(mode = "list", length = ncol(boot_samples))  
list_rsf <- vector(mode = "list", length = ncol(boot_samples))  
no_cr <- c()

# run the OOB validation
tic("Running the OOB validation")


for (i in 1:B) {
  
  cat("Started iteration i = ", i, "out of", B, "\n")
  # use the proper bootstrap column
  indices <- boot_samples[, i]
  oob <- eSTSmany$id[!eSTSmany$id %in% indices] # select oob patients
  
  #############################
  # training data
  #############################
  
  eSTSmany2 <- eSTSmany[indices, ]
  eSTSmany2$id <- 1:nrow(eSTSmany2)
  eSTSmany2$age <- scale(eSTSmany2$age)
  eSTSmany2$size <- scale(eSTSmany2$size)
  
  #############################
  # validation data
  #############################
  
  eSTSmany3 <- eSTSmany[oob, ]
  
  # if there is no competing cause on test data move to next i
  if (sum(eSTSmany3$status == 2) == 0)  {
    no_cr <- c(no_cr, i) # expanding vector with elements where OOB test set has no CR
    next 
  }
  
  eSTSmany3$age <- scale(eSTSmany3$age)
  eSTSmany3$size <- scale(eSTSmany3$size)
  
  #####################################
  # create the necessary long formats
  #####################################
  
  train_long_eSTSmany <- train_creator_eSTSmany(data = eSTSmany2)
  train_long_eSTSmany <- model.matrix(~., train_long_eSTSmany)[, -1] 
  
  train_x_eSTSmany <- train_long_eSTSmany[, c(4:17, 19)] # p = 14 variables + time interval
  dimnames(train_x_eSTSmany) <- NULL # the object must have empty dimnames
  train_y_eSTSmany <- train_long_eSTSmany[, (ncol(train_long_eSTSmany) - 2):ncol(train_long_eSTSmany)]
  dimnames(train_y_eSTSmany) <- NULL
  
  test_long_eSTSmany <- test_creator_eSTSmany(data = eSTSmany2)
  test_long_eSTSmany <- model.matrix(~., test_long_eSTSmany)[, -1]
  
  test_x_eSTSmany <- test_long_eSTSmany[, c(4:17, 19)] # p = 14 variables + time interval
  dimnames(test_x_eSTSmany) <- NULL # the object must have empty dimnames
  test_y_eSTSmany <- test_long_eSTSmany[, (ncol(train_long_eSTSmany) - 2):ncol(train_long_eSTSmany)]
  dimnames(test_y_eSTSmany) <- NULL
  
  valid_long_eSTSmany <- test_creator_eSTSmany(data = eSTSmany3)
  valid_long_eSTSmany <- model.matrix(~., valid_long_eSTSmany)[, -1]
  
  
  # build the model in nnet (validation for PLANNCR original)
  set.seed(12345)
  fit_nnet <- nnet(x = train_x_eSTSmany, y = train_y_eSTSmany,
                   size = 3, decay = 0.5, 
                   maxit = 500, trace = FALSE, softmax = TRUE)
  
  
  results1 <- measures_nnet_eSTSmany(trained_model = fit_nnet,
                                     datashort = eSTSmany3,
                                     datalong = valid_long_eSTSmany)
  list_nnet[[i]] <- results1
  
  
  # build the model in keras (validation for PLANNCR extended)
  # optimal parameters tuned with Brier score at 5 years
  # similarly the optimal paramters for AUC at 5 years can be selected 
  
  # set global default to never show metrics
  options(keras.view_metrics = FALSE)
  k_clear_session() # to avoid clutter from old models / layers in cross validation
  tensorflow::tf$random$set_seed(12345) # for update to TensorFlow 2.0
  # model with one hidden layer
  
  # Initialize a sequential model
  fit_keras1 <- keras_model_sequential()
  
  # Add layers to the model
  fit_keras1 %>%
    layer_dense(units = 12, activation = 'sigmoid',
                input_shape = ncol(train_x_eSTSmany)) %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 3, activation = 'softmax')
  
  # Compile the model
  fit_keras1 %>% compile(
    loss = 'categorical_crossentropy', # for multi-class classification problem
    optimizer = optimizer_sgd(learning_rate = 0.1, momentum = 0.8)
    #, metrics = c("accuracy")
  )
  early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 3)
  
  result_keras <- fit_keras1 %>% fit(
    train_x_eSTSmany,
    train_y_eSTSmany,
    epochs = 20,
    batch_size = 32,
    validation_data = list(test_x_eSTSmany, test_y_eSTSmany),
    class_weight = list("0" = 1, "1" = 1, "2" = 1),
    callbacks = c(early_stopping),
    verbose = 0 # don't display progress bar
  )
  
  results2 <- measures_keras_eSTSmany(trained_model = fit_keras1,
                                      datashort = eSTSmany3,
                                      datalong = valid_long_eSTSmany)
  
  list_keras[[i]] <- results2
  
  # 
  # # build Cause-Specific Cox
  
  cox.fit <- CSC(Hist(time,status)~ sex + age + size + margin + chemo
                 + grade + histology + depth + radioth, data=eSTSmany2)
  
  r1 <- predictRisk(cox.fit,newdata=as.data.table(eSTSmany3),times=0:11,cause=1)
  r2 <- predictRisk(cox.fit,newdata=as.data.table(eSTSmany3),times=0:11,cause=2)
  
  results3 <- measures_eSTSmany(risk_mat1 = r1, risk_mat2 = r2, datashort = eSTSmany3)
  list_csc[[i]] <- results3
  
  # build Fine - Gray
  fg.fit <- FGR(Hist(time,status)~ sex + age + size + margin + chemo
                + grade + histology + depth + radioth, data=eSTSmany2, cause=1)
  
  r1b <- predictRisk(fg.fit, newdata=as.data.table(eSTSmany3), times=0:11)
  
  fg.fit2 <- FGR(Hist(time,status)~ sex + age + size + margin + chemo
                 + grade + histology + depth + radioth, data=eSTSmany2, cause=2)
  r2b <- predictRisk(fg.fit2, newdata=as.data.table(eSTSmany3), times=0:11)
  
  results4 <- measures_eSTSmany(risk_mat1 = r1b, risk_mat2 = r2b, datashort = eSTSmany3)
  list_fg[[i]] <- results4
  
  # build RSFCR
  rsf.fit <- rfsrc(Surv(time, status) ~ sex + age + size + margin + chemo
                   + grade + histology + depth + radioth,
                   data = eSTSmany2, ntree = 1000, mtry = 5, nsplit = 2,
                   nodesize = 10, splitrule = "logrank", cause = 1, seed = 12345,
                   forest = TRUE, importance = FALSE)
  
  r1c <- predictRisk(rsf.fit, newdata=as.data.table(eSTSmany3), times=0:11, cause=1)
  r2c <- predictRisk(rsf.fit, newdata=as.data.table(eSTSmany3), times=0:11, cause=2)
  
  results5 <- measures_eSTSmany(risk_mat1 = r1c, risk_mat2 = r2c, datashort = eSTSmany3)
  list_rsf[[i]] <- results5
  
  
}

time_list <- toc()

cat((time_list$toc - time_list$tic) / 60, "minutes elapsed", "\n") 
# 34.54683 minutes elapsed

save(list_nnet, file = "results_nnet.RData")
save(list_keras, file = "results_keras_sigmoid_brier.RData")
save(list_rsf, file = "results_rsf_logrank.RData")
save(list_csc, file = "results_csc.RData")
save(list_fg, file = "results_fg.RData")



## end of the file



######################################################################
######################################################################
# file "bargraphs.R" 

# Used to to create bargraphs for all methods (for Brier score, AUC)
######################################################################
######################################################################

## start of the file


load("results_csc.RData")
load("results_fg.RData")
load("results_keras_sigmoid_brier.RData")
load("results_nnet.RData")
load("results_rsf_logrank.RData")
library(ggplot2)

# function that calculates the quantiles and the mean

summa <- function(x) {
  
  res1 <- round(quantile(x, probs = c(0, 0.025, 0.5, 0.975, 1), na.rm = TRUE), digits = 3)
  res2 <- round(mean(x, na.rm = TRUE), digits = 3)
  res3 <- round(sd(x, na.rm = TRUE), digits = 3)
  res <- c(res1[1:3], res2, res1[4:5], res3)
  names(res) <- c("Min.", "2.5% Qu.", "Median", "Mean", "97.5% Qu.", "Max.", "Sd.")
  return(res)
}



quant_low <- function(x) {
  res <- as.numeric(quantile(x, probs = 0.025, na.rm = TRUE))
  return(res)
}

quant_high <- function(x) {
  res <- as.numeric(quantile(x, probs = 0.975, na.rm = TRUE))
  return(res)
}




##########################################
# For cause-specific Cox model
##########################################


# cause 1
brier_csc_cause1 <- matrix(data = 0, nrow = length(list_csc), ncol = 11) # Brier score until 10 years
auc_csc_cause1 <- matrix(data = 0, nrow = length(list_csc), ncol = 11) # AUC until 10 years
calib_csc2y_cause1 <- vector(mode = "numeric", length = length(list_csc))
calib_csc5y_cause1 <- vector(mode = "numeric", length = length(list_csc))
calib_csc10y_cause1 <- vector(mode = "numeric", length = length(list_csc))

# cause 2
brier_csc_cause2 <- matrix(data = 0, nrow = length(list_csc), ncol = 11) # Brier score until 10 years
auc_csc_cause2 <- matrix(data = 0, nrow = length(list_csc), ncol = 11) # AUC until 10 years
calib_csc2y_cause2 <- vector(mode = "numeric", length = length(list_csc))
calib_csc5y_cause2 <- vector(mode = "numeric", length = length(list_csc))
calib_csc10y_cause2 <- vector(mode = "numeric", length = length(list_csc))


for (i in 1:length(list_csc)){
  
  
  brier_csc_cause1[i, ] <- list_csc[[i]]$Brier_cause1[1:11]
  brier_csc_cause2[i, ] <- list_csc[[i]]$Brier_cause2[1:11]
  auc_csc_cause1[i, ] <- list_csc[[i]]$AUC_cause1[1:11]
  auc_csc_cause2[i, ] <- list_csc[[i]]$AUC_cause2[1:11]
  
  calib_csc2y_cause1[i] <- list_csc[[i]]$MSE2_cause1
  calib_csc2y_cause2[i] <- list_csc[[i]]$MSE2_cause2
  calib_csc5y_cause1[i] <- list_csc[[i]]$MSE5_cause1
  calib_csc5y_cause2[i] <- list_csc[[i]]$MSE5_cause2
  calib_csc10y_cause1[i] <- list_csc[[i]]$MSE10_cause1
  calib_csc10y_cause2[i] <- list_csc[[i]]$MSE10_cause2 
  
}


##########################################
# For Fine-Gray model
##########################################


# cause 1
brier_fg_cause1 <- matrix(data = 0, nrow = length(list_fg), ncol = 11) # Brier score until 10 years
auc_fg_cause1 <- matrix(data = 0, nrow = length(list_fg), ncol = 11) # AUC until 10 years
calib_fg2y_cause1 <- vector(mode = "numeric", length = length(list_fg))
calib_fg5y_cause1 <- vector(mode = "numeric", length = length(list_fg))
calib_fg10y_cause1 <- vector(mode = "numeric", length = length(list_fg))

# cause 2
brier_fg_cause2 <- matrix(data = 0, nrow = length(list_fg), ncol = 11) # Brier score until 10 years
auc_fg_cause2 <- matrix(data = 0, nrow = length(list_fg), ncol = 11) # AUC until 10 years
calib_fg2y_cause2 <- vector(mode = "numeric", length = length(list_fg))
calib_fg5y_cause2 <- vector(mode = "numeric", length = length(list_fg))
calib_fg10y_cause2 <- vector(mode = "numeric", length = length(list_fg))



for (i in 1:length(list_fg)){
  
  
  brier_fg_cause1[i, ] <- list_fg[[i]]$Brier_cause1[1:11]
  brier_fg_cause2[i, ] <- list_fg[[i]]$Brier_cause2[1:11]
  auc_fg_cause1[i, ] <- list_fg[[i]]$AUC_cause1[1:11]
  auc_fg_cause2[i, ] <- list_fg[[i]]$AUC_cause2[1:11]
  
  calib_fg2y_cause1[i] <- list_fg[[i]]$MSE2_cause1
  calib_fg2y_cause2[i] <- list_fg[[i]]$MSE2_cause2
  calib_fg5y_cause1[i] <- list_fg[[i]]$MSE5_cause1
  calib_fg5y_cause2[i] <- list_fg[[i]]$MSE5_cause2
  calib_fg10y_cause1[i] <- list_fg[[i]]$MSE10_cause1
  calib_fg10y_cause2[i] <- list_fg[[i]]$MSE10_cause2 
  
}




##########################################
# For PLANNCR original
##########################################

# cause 1
brier_nnet_cause1 <- matrix(data = 0, nrow = length(list_nnet), ncol = 11) # Brier score until 10 years
auc_nnet_cause1 <- matrix(data = 0, nrow = length(list_nnet), ncol = 11) # AUC until 10 years
calib_nnet2y_cause1 <- vector(mode = "numeric", length = length(list_nnet))
calib_nnet5y_cause1 <- vector(mode = "numeric", length = length(list_nnet))
calib_nnet10y_cause1 <- vector(mode = "numeric", length = length(list_nnet))

# cause 2
brier_nnet_cause2 <- matrix(data = 0, nrow = length(list_nnet), ncol = 11) # Brier score until 10 years
auc_nnet_cause2 <- matrix(data = 0, nrow = length(list_nnet), ncol = 11) # AUC until 10 years
calib_nnet2y_cause2 <- vector(mode = "numeric", length = length(list_nnet))
calib_nnet5y_cause2 <- vector(mode = "numeric", length = length(list_nnet))
calib_nnet10y_cause2 <- vector(mode = "numeric", length = length(list_nnet))


for (i in 1:length(list_nnet)){
  
  
  brier_nnet_cause1[i, ] <- list_nnet[[i]]$Brier_cause1[1:11]
  brier_nnet_cause2[i, ] <- list_nnet[[i]]$Brier_cause2[1:11]
  auc_nnet_cause1[i, ] <- list_nnet[[i]]$AUC_cause1[1:11]
  auc_nnet_cause2[i, ] <- list_nnet[[i]]$AUC_cause2[1:11]
  
  calib_nnet2y_cause1[i] <- list_nnet[[i]]$MSE2_cause1
  calib_nnet2y_cause2[i] <- list_nnet[[i]]$MSE2_cause2
  calib_nnet5y_cause1[i] <- list_nnet[[i]]$MSE5_cause1
  calib_nnet5y_cause2[i] <- list_nnet[[i]]$MSE5_cause2
  calib_nnet10y_cause1[i] <- list_nnet[[i]]$MSE10_cause1
  calib_nnet10y_cause2[i] <- list_nnet[[i]]$MSE10_cause2 
  
}

##########################################
# For PLANNCR extended
##########################################

# cause 1
brier_keras_cause1 <- matrix(data = 0, nrow = length(list_keras), ncol = 11) # Brier score until 10 years
auc_keras_cause1 <- matrix(data = 0, nrow = length(list_keras), ncol = 11) # AUC until 10 years
calib_keras2y_cause1 <- vector(mode = "numeric", length = length(list_keras))
calib_keras5y_cause1 <- vector(mode = "numeric", length = length(list_keras))
calib_keras10y_cause1 <- vector(mode = "numeric", length = length(list_keras))

# cause 2
brier_keras_cause2 <- matrix(data = 0, nrow = length(list_keras), ncol = 11) # Brier score until 10 years
auc_keras_cause2 <- matrix(data = 0, nrow = length(list_keras), ncol = 11) # AUC until 10 years
calib_keras2y_cause2 <- vector(mode = "numeric", length = length(list_keras))
calib_keras5y_cause2 <- vector(mode = "numeric", length = length(list_keras))
calib_keras10y_cause2 <- vector(mode = "numeric", length = length(list_keras))


for (i in 1:length(list_keras)){
  
  
  brier_keras_cause1[i, ] <- list_keras[[i]]$Brier_cause1[1:11]
  brier_keras_cause2[i, ] <- list_keras[[i]]$Brier_cause2[1:11]
  auc_keras_cause1[i, ] <- list_keras[[i]]$AUC_cause1[1:11]
  auc_keras_cause2[i, ] <- list_keras[[i]]$AUC_cause2[1:11]
  
  calib_keras2y_cause1[i] <- list_keras[[i]]$MSE2_cause1
  calib_keras2y_cause2[i] <- list_keras[[i]]$MSE2_cause2
  calib_keras5y_cause1[i] <- list_keras[[i]]$MSE5_cause1
  calib_keras5y_cause2[i] <- list_keras[[i]]$MSE5_cause2
  calib_keras10y_cause1[i] <- list_keras[[i]]$MSE10_cause1
  calib_keras10y_cause2[i] <- list_keras[[i]]$MSE10_cause2 
  
}


##########################################
# For RSFCR
##########################################

# cause 1
brier_rsf_cause1 <- matrix(data = 0, nrow = length(list_rsf), ncol = 11) # Brier score until 10 years
auc_rsf_cause1 <- matrix(data = 0, nrow = length(list_rsf), ncol = 11) # AUC until 10 years
calib_rsf2y_cause1 <- vector(mode = "numeric", length = length(list_rsf))
calib_rsf5y_cause1 <- vector(mode = "numeric", length = length(list_rsf))
calib_rsf10y_cause1 <- vector(mode = "numeric", length = length(list_rsf))

# cause 2
brier_rsf_cause2 <- matrix(data = 0, nrow = length(list_rsf), ncol = 11) # Brier score until 10 years
auc_rsf_cause2 <- matrix(data = 0, nrow = length(list_rsf), ncol = 11) # AUC until 10 years
calib_rsf2y_cause2 <- vector(mode = "numeric", length = length(list_rsf))
calib_rsf5y_cause2 <- vector(mode = "numeric", length = length(list_rsf))
calib_rsf10y_cause2 <- vector(mode = "numeric", length = length(list_rsf))



for (i in 1:length(list_rsf)){
  
  brier_rsf_cause1[i, ] <- list_rsf[[i]]$Brier_cause1[1:11]
  brier_rsf_cause2[i, ] <- list_rsf[[i]]$Brier_cause2[1:11]
  auc_rsf_cause1[i, ] <- list_rsf[[i]]$AUC_cause1[1:11]
  auc_rsf_cause2[i, ] <- list_rsf[[i]]$AUC_cause2[1:11]
  
  calib_rsf2y_cause1[i] <- list_rsf[[i]]$MSE2_cause1
  calib_rsf2y_cause2[i] <- list_rsf[[i]]$MSE2_cause2
  calib_rsf5y_cause1[i] <- list_rsf[[i]]$MSE5_cause1
  calib_rsf5y_cause2[i] <- list_rsf[[i]]$MSE5_cause2
  calib_rsf10y_cause1[i] <- list_rsf[[i]]$MSE10_cause1
  calib_rsf10y_cause2[i] <- list_rsf[[i]]$MSE10_cause2 
  
}


# Bargraph for Brier score - AUC (for the event of interest: disease progression)

# at 2 years

line1 <- data.frame(brier = mean(brier_csc_cause1[, 3], na.rm = TRUE),
                    brier_left = quant_low(brier_csc_cause1[, 3]),
                    brier_right = quant_high(brier_csc_cause1[, 3]),
                    auc = mean(auc_csc_cause1[, 3], na.rm = TRUE),
                    auc_left = quant_low(auc_csc_cause1[, 3]),
                    auc_right = quant_high(auc_csc_cause1[, 3]),
                    Model = "Cause-specific Cox",
                    Time = "2 years")


line2 <- data.frame(brier = mean(brier_fg_cause1[, 3], na.rm = TRUE),
                    brier_left = quant_low(brier_fg_cause1[, 3]),
                    brier_right = quant_high(brier_fg_cause1[, 3]),
                    auc = mean(auc_fg_cause1[, 3], na.rm = TRUE),
                    auc_left = quant_low(auc_fg_cause1[, 3]),
                    auc_right = quant_high(auc_fg_cause1[, 3]),
                    Model = "Fine-Gray",
                    Time = "2 years")


line3 <- data.frame(brier = mean(brier_nnet_cause1[, 3], na.rm = TRUE),
                    brier_left = quant_low(brier_nnet_cause1[, 3]),
                    brier_right = quant_high(brier_nnet_cause1[, 3]),
                    auc = mean(auc_nnet_cause1[, 3], na.rm = TRUE),
                    auc_left = quant_low(auc_nnet_cause1[, 3]),
                    auc_right = quant_high(auc_nnet_cause1[, 3]),
                    Model = "PLANNCR original",
                    Time = "2 years")

line4 <- data.frame(brier = mean(brier_keras_cause1[, 3], na.rm = TRUE),
                    brier_left = quant_low(brier_keras_cause1[, 3]),
                    brier_right = quant_high(brier_keras_cause1[, 3]),
                    auc = mean(auc_keras_cause1[, 3], na.rm = TRUE),
                    auc_left = quant_low(auc_keras_cause1[, 3]),
                    auc_right = quant_high(auc_keras_cause1[, 3]),
                    Model = "PLANNCR extended",
                    Time = "2 years")

line5 <- data.frame(brier = mean(brier_rsf_cause1[, 3], na.rm = TRUE),
                    brier_left = quant_low(brier_rsf_cause1[, 3]),
                    brier_right = quant_high(brier_rsf_cause1[, 3]),
                    auc = mean(auc_rsf_cause1[, 3], na.rm = TRUE),
                    auc_left = quant_low(auc_rsf_cause1[, 3]),
                    auc_right = quant_high(auc_rsf_cause1[, 3]),
                    Model = "RSFCR",
                    Time = "2 years")



# at 5 years

line6 <- data.frame(brier = mean(brier_csc_cause1[, 6], na.rm = TRUE),
                    brier_left = quant_low(brier_csc_cause1[, 6]),
                    brier_right = quant_high(brier_csc_cause1[, 6]),
                    auc = mean(auc_csc_cause1[, 6], na.rm = TRUE),
                    auc_left = quant_low(auc_csc_cause1[, 6]),
                    auc_right = quant_high(auc_csc_cause1[, 6]),
                    Model = "Cause-specific Cox",
                    Time = "5 years")


line7 <- data.frame(brier = mean(brier_fg_cause1[, 6], na.rm = TRUE),
                    brier_left = quant_low(brier_fg_cause1[, 6]),
                    brier_right = quant_high(brier_fg_cause1[, 6]),
                    auc = mean(auc_fg_cause1[, 6], na.rm = TRUE),
                    auc_left = quant_low(auc_fg_cause1[, 6]),
                    auc_right = quant_high(auc_fg_cause1[, 6]),
                    Model = "Fine-Gray",
                    Time = "5 years")


line8 <- data.frame(brier = mean(brier_nnet_cause1[, 6], na.rm = TRUE),
                    brier_left = quant_low(brier_nnet_cause1[, 6]),
                    brier_right = quant_high(brier_nnet_cause1[, 6]),
                    auc = mean(auc_nnet_cause1[, 6], na.rm = TRUE),
                    auc_left = quant_low(auc_nnet_cause1[, 6]),
                    auc_right = quant_high(auc_nnet_cause1[, 6]),
                    Model = "PLANNCR original",
                    Time = "5 years")

line9 <- data.frame(brier = mean(brier_keras_cause1[, 6], na.rm = TRUE),
                    brier_left = quant_low(brier_keras_cause1[, 6]),
                    brier_right = quant_high(brier_keras_cause1[, 6]),
                    auc = mean(auc_keras_cause1[, 6], na.rm = TRUE),
                    auc_left = quant_low(auc_keras_cause1[, 6]),
                    auc_right = quant_high(auc_keras_cause1[, 6]),
                    Model = "PLANNCR extended",
                    Time = "5 years")

line10 <- data.frame(brier = mean(brier_rsf_cause1[, 6], na.rm = TRUE),
                     brier_left = quant_low(brier_rsf_cause1[, 6]),
                     brier_right = quant_high(brier_rsf_cause1[, 6]),
                     auc = mean(auc_rsf_cause1[, 6], na.rm = TRUE),
                     auc_left = quant_low(auc_rsf_cause1[, 6]),
                     auc_right = quant_high(auc_rsf_cause1[, 6]),
                     Model = "RSFCR",
                     Time = "5 years")



# at 10 years

line11 <- data.frame(brier = mean(brier_csc_cause1[, 11], na.rm = TRUE),
                     brier_left = quant_low(brier_csc_cause1[, 11]),
                     brier_right = quant_high(brier_csc_cause1[, 11]),
                     auc = mean(auc_csc_cause1[, 11], na.rm = TRUE),
                     auc_left = quant_low(auc_csc_cause1[, 11]),
                     auc_right = quant_high(auc_csc_cause1[, 11]),
                     Model = "Cause-specific Cox",
                     Time = "10 years")


line12 <- data.frame(brier = mean(brier_fg_cause1[, 11], na.rm = TRUE),
                     brier_left = quant_low(brier_fg_cause1[, 11]),
                     brier_right = quant_high(brier_fg_cause1[, 11]),
                     auc = mean(auc_fg_cause1[, 11], na.rm = TRUE),
                     auc_left = quant_low(auc_fg_cause1[, 11]),
                     auc_right = quant_high(auc_fg_cause1[, 11]),
                     Model = "Fine-Gray",
                     Time = "10 years")


line13 <- data.frame(brier = mean(brier_nnet_cause1[, 11], na.rm = TRUE),
                     brier_left = quant_low(brier_nnet_cause1[, 11]),
                     brier_right = quant_high(brier_nnet_cause1[, 11]),
                     auc = mean(auc_nnet_cause1[, 11], na.rm = TRUE),
                     auc_left = quant_low(auc_nnet_cause1[, 11]),
                     auc_right = quant_high(auc_nnet_cause1[, 11]),
                     Model = "PLANNCR original",
                     Time = "10 years")

line14 <- data.frame(brier = mean(brier_keras_cause1[, 11], na.rm = TRUE),
                     brier_left = quant_low(brier_keras_cause1[, 11]),
                     brier_right = quant_high(brier_keras_cause1[, 11]),
                     auc = mean(auc_keras_cause1[, 11], na.rm = TRUE),
                     auc_left = quant_low(auc_keras_cause1[, 11]),
                     auc_right = quant_high(auc_keras_cause1[, 11]),
                     Model = "PLANNCR extended",
                     Time = "10 years")

line15 <- data.frame(brier = mean(brier_rsf_cause1[, 11], na.rm = TRUE),
                     brier_left = quant_low(brier_rsf_cause1[, 11]),
                     brier_right = quant_high(brier_rsf_cause1[, 11]),
                     auc = mean(auc_rsf_cause1[, 11], na.rm = TRUE),
                     auc_left = quant_low(auc_rsf_cause1[, 11]),
                     auc_right = quant_high(auc_rsf_cause1[, 11]),
                     Model = "RSFCR",
                     Time = "10 years")




df <- rbind(line1, line2, line3, line4, line5, line6, line7, line8,
            line9, line10, line11, line12, line13, line14, line15)

df$Model <- as.factor(df$Model)
df$Model <- factor(df$Model, levels = c("Cause-specific Cox", "Fine-Gray", "PLANNCR original",
                                        "PLANNCR extended", "RSFCR"))
df$Time <- as.factor(df$Time)
df$Time <- factor(df$Time, levels = c("2 years", "5 years", "10 years"))



# For the Brier scores

# Use 95% confidence intervals 
brier_c1 <- ggplot(df, aes(x=Time, y=brier, fill = Model)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  labs(title = " ", x = "Time since surgery",
       y = "Brier score (for disease progression)", fill = "Model") +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25), limits = c(0, 0.26)) +
  geom_errorbar(aes(ymin=brier_left, ymax=brier_right),
                size = 0.5,    # Thinner lines
                width = 0.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top") # Remove legend by setting to "none" 
# + theme_classic() 


# For the AUC

auc_c1 <- ggplot(df, aes(x=Time, y=auc, fill = Model)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  labs(title = " ", x = "Time since surgery",
       y = "AUC (for disease progression)", fill = "Model") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1)) +
  geom_errorbar(aes(ymin=auc_left, ymax=auc_right),
                size = 0.5,    # Thinner lines
                width = 0.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top") # Remove legend by setting to "none" 
# + theme_classic() 


gridExtra::grid.arrange(brier_c1, auc_c1, ncol = 2)



# for supplementary material


# Bargraph for Brier score - AUC (for the competing event: death)

# at 2 years

line16 <- data.frame(brier = mean(brier_csc_cause2[, 3], na.rm = TRUE),
                     brier_left = quant_low(brier_csc_cause2[, 3]),
                     brier_right = quant_high(brier_csc_cause2[, 3]),
                     auc = mean(auc_csc_cause2[, 3], na.rm = TRUE),
                     auc_left = quant_low(auc_csc_cause2[, 3]),
                     auc_right = quant_high(auc_csc_cause2[, 3]),
                     Model = "Cause-specific Cox",
                     Time = "2 years")


line17 <- data.frame(brier = mean(brier_fg_cause2[, 3], na.rm = TRUE),
                     brier_left = quant_low(brier_fg_cause2[, 3]),
                     brier_right = quant_high(brier_fg_cause2[, 3]),
                     auc = mean(auc_fg_cause2[, 3], na.rm = TRUE),
                     auc_left = quant_low(auc_fg_cause2[, 3]),
                     auc_right = quant_high(auc_fg_cause2[, 3]),
                     Model = "Fine-Gray",
                     Time = "2 years")


line18 <- data.frame(brier = mean(brier_nnet_cause2[, 3], na.rm = TRUE),
                     brier_left = quant_low(brier_nnet_cause2[, 3]),
                     brier_right = quant_high(brier_nnet_cause2[, 3]),
                     auc = mean(auc_nnet_cause2[, 3], na.rm = TRUE),
                     auc_left = quant_low(auc_nnet_cause2[, 3]),
                     auc_right = quant_high(auc_nnet_cause2[, 3]),
                     Model = "PLANNCR original",
                     Time = "2 years")

line19 <- data.frame(brier = mean(brier_keras_cause2[, 3], na.rm = TRUE),
                     brier_left = quant_low(brier_keras_cause2[, 3]),
                     brier_right = quant_high(brier_keras_cause2[, 3]),
                     auc = mean(auc_keras_cause2[, 3], na.rm = TRUE),
                     auc_left = quant_low(auc_keras_cause2[, 3]),
                     auc_right = quant_high(auc_keras_cause2[, 3]),
                     Model = "PLANNCR extended",
                     Time = "2 years")

line20 <- data.frame(brier = mean(brier_rsf_cause2[, 3], na.rm = TRUE),
                     brier_left = quant_low(brier_rsf_cause2[, 3]),
                     brier_right = quant_high(brier_rsf_cause2[, 3]),
                     auc = mean(auc_rsf_cause2[, 3], na.rm = TRUE),
                     auc_left = quant_low(auc_rsf_cause2[, 3]),
                     auc_right = quant_high(auc_rsf_cause2[, 3]),
                     Model = "RSFCR",
                     Time = "2 years")



# at 5 years

line21 <- data.frame(brier = mean(brier_csc_cause2[, 6], na.rm = TRUE),
                     brier_left = quant_low(brier_csc_cause2[, 6]),
                     brier_right = quant_high(brier_csc_cause2[, 6]),
                     auc = mean(auc_csc_cause2[, 6], na.rm = TRUE),
                     auc_left = quant_low(auc_csc_cause2[, 6]),
                     auc_right = quant_high(auc_csc_cause2[, 6]),
                     Model = "Cause-specific Cox",
                     Time = "5 years")


line22 <- data.frame(brier = mean(brier_fg_cause2[, 6], na.rm = TRUE),
                     brier_left = quant_low(brier_fg_cause2[, 6]),
                     brier_right = quant_high(brier_fg_cause2[, 6]),
                     auc = mean(auc_fg_cause2[, 6], na.rm = TRUE),
                     auc_left = quant_low(auc_fg_cause2[, 6]),
                     auc_right = quant_high(auc_fg_cause2[, 6]),
                     Model = "Fine-Gray",
                     Time = "5 years")


line23 <- data.frame(brier = mean(brier_nnet_cause2[, 6], na.rm = TRUE),
                     brier_left = quant_low(brier_nnet_cause2[, 6]),
                     brier_right = quant_high(brier_nnet_cause2[, 6]),
                     auc = mean(auc_nnet_cause2[, 6], na.rm = TRUE),
                     auc_left = quant_low(auc_nnet_cause2[, 6]),
                     auc_right = quant_high(auc_nnet_cause2[, 6]),
                     Model = "PLANNCR original",
                     Time = "5 years")

line24 <- data.frame(brier = mean(brier_keras_cause2[, 6], na.rm = TRUE),
                     brier_left = quant_low(brier_keras_cause2[, 6]),
                     brier_right = quant_high(brier_keras_cause2[, 6]),
                     auc = mean(auc_keras_cause2[, 6], na.rm = TRUE),
                     auc_left = quant_low(auc_keras_cause2[, 6]),
                     auc_right = quant_high(auc_keras_cause2[, 6]),
                     Model = "PLANNCR extended",
                     Time = "5 years")

line25 <- data.frame(brier = mean(brier_rsf_cause2[, 6], na.rm = TRUE),
                     brier_left = quant_low(brier_rsf_cause2[, 6]),
                     brier_right = quant_high(brier_rsf_cause2[, 6]),
                     auc = mean(auc_rsf_cause2[, 6], na.rm = TRUE),
                     auc_left = quant_low(auc_rsf_cause2[, 6]),
                     auc_right = quant_high(auc_rsf_cause2[, 6]),
                     Model = "RSFCR",
                     Time = "5 years")



# at 10 years

line26 <- data.frame(brier = mean(brier_csc_cause2[, 11], na.rm = TRUE),
                     brier_left = quant_low(brier_csc_cause2[, 11]),
                     brier_right = quant_high(brier_csc_cause2[, 11]),
                     auc = mean(auc_csc_cause2[, 11], na.rm = TRUE),
                     auc_left = quant_low(auc_csc_cause2[, 11]),
                     auc_right = quant_high(auc_csc_cause2[, 11]),
                     Model = "Cause-specific Cox",
                     Time = "10 years")


line27 <- data.frame(brier = mean(brier_fg_cause2[, 11], na.rm = TRUE),
                     brier_left = quant_low(brier_fg_cause2[, 11]),
                     brier_right = quant_high(brier_fg_cause2[, 11]),
                     auc = mean(auc_fg_cause2[, 11], na.rm = TRUE),
                     auc_left = quant_low(auc_fg_cause2[, 11]),
                     auc_right = quant_high(auc_fg_cause2[, 11]),
                     Model = "Fine-Gray",
                     Time = "10 years")


line28 <- data.frame(brier = mean(brier_nnet_cause2[, 11], na.rm = TRUE),
                     brier_left = quant_low(brier_nnet_cause2[, 11]),
                     brier_right = quant_high(brier_nnet_cause2[, 11]),
                     auc = mean(auc_nnet_cause2[, 11], na.rm = TRUE),
                     auc_left = quant_low(auc_nnet_cause2[, 11]),
                     auc_right = quant_high(auc_nnet_cause2[, 11]),
                     Model = "PLANNCR original",
                     Time = "10 years")

line29 <- data.frame(brier = mean(brier_keras_cause2[, 11], na.rm = TRUE),
                     brier_left = quant_low(brier_keras_cause2[, 11]),
                     brier_right = quant_high(brier_keras_cause2[, 11]),
                     auc = mean(auc_keras_cause2[, 11], na.rm = TRUE),
                     auc_left = quant_low(auc_keras_cause2[, 11]),
                     auc_right = quant_high(auc_keras_cause2[, 11]),
                     Model = "PLANNCR extended",
                     Time = "10 years")

line30 <- data.frame(brier = mean(brier_rsf_cause2[, 11], na.rm = TRUE),
                     brier_left = quant_low(brier_rsf_cause2[, 11]),
                     brier_right = quant_high(brier_rsf_cause2[, 11]),
                     auc = mean(auc_rsf_cause2[, 11], na.rm = TRUE),
                     auc_left = quant_low(auc_rsf_cause2[, 11]),
                     auc_right = quant_high(auc_rsf_cause2[, 11]),
                     Model = "RSFCR",
                     Time = "10 years")




df2 <- rbind(line16, line17, line18, line19, line20, line21, line22, line23,
             line24, line25, line26, line27, line28, line29, line30)

df2$Model <- as.factor(df2$Model)
df2$Model <- factor(df2$Model, levels = c("Cause-specific Cox", "Fine-Gray", "PLANNCR original",
                                          "PLANNCR extended", "RSFCR"))
df2$Time <- as.factor(df2$Time)
df2$Time <- factor(df2$Time, levels = c("2 years", "5 years", "10 years"))



# For the Brier scores

# Use 95% confidence intervals 
brier_c2 <- ggplot(df2, aes(x=Time, y=brier, fill = Model)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  labs(title = " ", x = "Time since surgery",
       y = "Brier score (for death)", fill = "Model") +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25), limits = c(0, 0.25)) +
  geom_errorbar(aes(ymin=brier_left, ymax=brier_right),
                size = 0.5,    # Thinner lines
                width = 0.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top") # Remove legend by setting to "none" 
# + theme_classic() 


# For the AUC

auc_c2 <- ggplot(df, aes(x=Time, y=auc, fill = Model)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  labs(title = " ", x = "Time since surgery",
       y = "AUC (for death)", fill = "Model") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1)) +
  geom_errorbar(aes(ymin=auc_left, ymax=auc_right),
                size = 0.5,    # Thinner lines
                width = 0.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top") # Remove legend by setting to "none" 
# + theme_classic() 


gridExtra::grid.arrange(brier_c2, auc_c2, ncol = 2)


## end of the file



######################################################################
######################################################################
# file "boxplots.R" 

# Used to to create bargraphs for all methods (for miscalibration)
######################################################################
######################################################################

## start of the file

load("results_csc.RData")
load("results_fg.RData")
load("results_keras_sigmoid_brier.RData")
load("results_nnet.RData")
load("results_rsf_logrank.RData")
library(ggplot2)

##########################################
# For cause-specific Cox model
##########################################


# cause 1
brier_csc_cause1 <- matrix(data = 0, nrow = length(list_csc), ncol = 11) # Brier score until 10 years
auc_csc_cause1 <- matrix(data = 0, nrow = length(list_csc), ncol = 11) # AUC until 10 years
calib_csc2y_cause1 <- vector(mode = "numeric", length = length(list_csc))
calib_csc5y_cause1 <- vector(mode = "numeric", length = length(list_csc))
calib_csc10y_cause1 <- vector(mode = "numeric", length = length(list_csc))

# cause 2
brier_csc_cause2 <- matrix(data = 0, nrow = length(list_csc), ncol = 11) # Brier score until 10 years
auc_csc_cause2 <- matrix(data = 0, nrow = length(list_csc), ncol = 11) # AUC until 10 years
calib_csc2y_cause2 <- vector(mode = "numeric", length = length(list_csc))
calib_csc5y_cause2 <- vector(mode = "numeric", length = length(list_csc))
calib_csc10y_cause2 <- vector(mode = "numeric", length = length(list_csc))


for (i in 1:length(list_csc)){
  
  
  brier_csc_cause1[i, ] <- list_csc[[i]]$Brier_cause1[1:11]
  brier_csc_cause2[i, ] <- list_csc[[i]]$Brier_cause2[1:11]
  auc_csc_cause1[i, ] <- list_csc[[i]]$AUC_cause1[1:11]
  auc_csc_cause2[i, ] <- list_csc[[i]]$AUC_cause2[1:11]
  
  calib_csc2y_cause1[i] <- list_csc[[i]]$MSE2_cause1
  calib_csc2y_cause2[i] <- list_csc[[i]]$MSE2_cause2
  calib_csc5y_cause1[i] <- list_csc[[i]]$MSE5_cause1
  calib_csc5y_cause2[i] <- list_csc[[i]]$MSE5_cause2
  calib_csc10y_cause1[i] <- list_csc[[i]]$MSE10_cause1
  calib_csc10y_cause2[i] <- list_csc[[i]]$MSE10_cause2 
  
}


##########################################
# For Fine-Gray model
##########################################


# cause 1
brier_fg_cause1 <- matrix(data = 0, nrow = length(list_fg), ncol = 11) # Brier score until 10 years
auc_fg_cause1 <- matrix(data = 0, nrow = length(list_fg), ncol = 11) # AUC until 10 years
calib_fg2y_cause1 <- vector(mode = "numeric", length = length(list_fg))
calib_fg5y_cause1 <- vector(mode = "numeric", length = length(list_fg))
calib_fg10y_cause1 <- vector(mode = "numeric", length = length(list_fg))

# cause 2
brier_fg_cause2 <- matrix(data = 0, nrow = length(list_fg), ncol = 11) # Brier score until 10 years
auc_fg_cause2 <- matrix(data = 0, nrow = length(list_fg), ncol = 11) # AUC until 10 years
calib_fg2y_cause2 <- vector(mode = "numeric", length = length(list_fg))
calib_fg5y_cause2 <- vector(mode = "numeric", length = length(list_fg))
calib_fg10y_cause2 <- vector(mode = "numeric", length = length(list_fg))



for (i in 1:length(list_fg)){
  
  
  brier_fg_cause1[i, ] <- list_fg[[i]]$Brier_cause1[1:11]
  brier_fg_cause2[i, ] <- list_fg[[i]]$Brier_cause2[1:11]
  auc_fg_cause1[i, ] <- list_fg[[i]]$AUC_cause1[1:11]
  auc_fg_cause2[i, ] <- list_fg[[i]]$AUC_cause2[1:11]
  
  calib_fg2y_cause1[i] <- list_fg[[i]]$MSE2_cause1
  calib_fg2y_cause2[i] <- list_fg[[i]]$MSE2_cause2
  calib_fg5y_cause1[i] <- list_fg[[i]]$MSE5_cause1
  calib_fg5y_cause2[i] <- list_fg[[i]]$MSE5_cause2
  calib_fg10y_cause1[i] <- list_fg[[i]]$MSE10_cause1
  calib_fg10y_cause2[i] <- list_fg[[i]]$MSE10_cause2 
  
}




##########################################
# For PLANNCR original
##########################################

# cause 1
brier_nnet_cause1 <- matrix(data = 0, nrow = length(list_nnet), ncol = 11) # Brier score until 10 years
auc_nnet_cause1 <- matrix(data = 0, nrow = length(list_nnet), ncol = 11) # AUC until 10 years
calib_nnet2y_cause1 <- vector(mode = "numeric", length = length(list_nnet))
calib_nnet5y_cause1 <- vector(mode = "numeric", length = length(list_nnet))
calib_nnet10y_cause1 <- vector(mode = "numeric", length = length(list_nnet))

# cause 2
brier_nnet_cause2 <- matrix(data = 0, nrow = length(list_nnet), ncol = 11) # Brier score until 10 years
auc_nnet_cause2 <- matrix(data = 0, nrow = length(list_nnet), ncol = 11) # AUC until 10 years
calib_nnet2y_cause2 <- vector(mode = "numeric", length = length(list_nnet))
calib_nnet5y_cause2 <- vector(mode = "numeric", length = length(list_nnet))
calib_nnet10y_cause2 <- vector(mode = "numeric", length = length(list_nnet))


for (i in 1:length(list_nnet)){
  
  
  brier_nnet_cause1[i, ] <- list_nnet[[i]]$Brier_cause1[1:11]
  brier_nnet_cause2[i, ] <- list_nnet[[i]]$Brier_cause2[1:11]
  auc_nnet_cause1[i, ] <- list_nnet[[i]]$AUC_cause1[1:11]
  auc_nnet_cause2[i, ] <- list_nnet[[i]]$AUC_cause2[1:11]
  
  calib_nnet2y_cause1[i] <- list_nnet[[i]]$MSE2_cause1
  calib_nnet2y_cause2[i] <- list_nnet[[i]]$MSE2_cause2
  calib_nnet5y_cause1[i] <- list_nnet[[i]]$MSE5_cause1
  calib_nnet5y_cause2[i] <- list_nnet[[i]]$MSE5_cause2
  calib_nnet10y_cause1[i] <- list_nnet[[i]]$MSE10_cause1
  calib_nnet10y_cause2[i] <- list_nnet[[i]]$MSE10_cause2 
  
}

##########################################
# For PLANNCR extended
##########################################

# cause 1
brier_keras_cause1 <- matrix(data = 0, nrow = length(list_keras), ncol = 11) # Brier score until 10 years
auc_keras_cause1 <- matrix(data = 0, nrow = length(list_keras), ncol = 11) # AUC until 10 years
calib_keras2y_cause1 <- vector(mode = "numeric", length = length(list_keras))
calib_keras5y_cause1 <- vector(mode = "numeric", length = length(list_keras))
calib_keras10y_cause1 <- vector(mode = "numeric", length = length(list_keras))

# cause 2
brier_keras_cause2 <- matrix(data = 0, nrow = length(list_keras), ncol = 11) # Brier score until 10 years
auc_keras_cause2 <- matrix(data = 0, nrow = length(list_keras), ncol = 11) # AUC until 10 years
calib_keras2y_cause2 <- vector(mode = "numeric", length = length(list_keras))
calib_keras5y_cause2 <- vector(mode = "numeric", length = length(list_keras))
calib_keras10y_cause2 <- vector(mode = "numeric", length = length(list_keras))


for (i in 1:length(list_keras)){
  
  
  brier_keras_cause1[i, ] <- list_keras[[i]]$Brier_cause1[1:11]
  brier_keras_cause2[i, ] <- list_keras[[i]]$Brier_cause2[1:11]
  auc_keras_cause1[i, ] <- list_keras[[i]]$AUC_cause1[1:11]
  auc_keras_cause2[i, ] <- list_keras[[i]]$AUC_cause2[1:11]
  
  calib_keras2y_cause1[i] <- list_keras[[i]]$MSE2_cause1
  calib_keras2y_cause2[i] <- list_keras[[i]]$MSE2_cause2
  calib_keras5y_cause1[i] <- list_keras[[i]]$MSE5_cause1
  calib_keras5y_cause2[i] <- list_keras[[i]]$MSE5_cause2
  calib_keras10y_cause1[i] <- list_keras[[i]]$MSE10_cause1
  calib_keras10y_cause2[i] <- list_keras[[i]]$MSE10_cause2 
  
}


##########################################
# For RSFCR
##########################################

# cause 1
brier_rsf_cause1 <- matrix(data = 0, nrow = length(list_rsf), ncol = 11) # Brier score until 10 years
auc_rsf_cause1 <- matrix(data = 0, nrow = length(list_rsf), ncol = 11) # AUC until 10 years
calib_rsf2y_cause1 <- vector(mode = "numeric", length = length(list_rsf))
calib_rsf5y_cause1 <- vector(mode = "numeric", length = length(list_rsf))
calib_rsf10y_cause1 <- vector(mode = "numeric", length = length(list_rsf))

# cause 2
brier_rsf_cause2 <- matrix(data = 0, nrow = length(list_rsf), ncol = 11) # Brier score until 10 years
auc_rsf_cause2 <- matrix(data = 0, nrow = length(list_rsf), ncol = 11) # AUC until 10 years
calib_rsf2y_cause2 <- vector(mode = "numeric", length = length(list_rsf))
calib_rsf5y_cause2 <- vector(mode = "numeric", length = length(list_rsf))
calib_rsf10y_cause2 <- vector(mode = "numeric", length = length(list_rsf))



for (i in 1:length(list_rsf)){
  
  brier_rsf_cause1[i, ] <- list_rsf[[i]]$Brier_cause1[1:11]
  brier_rsf_cause2[i, ] <- list_rsf[[i]]$Brier_cause2[1:11]
  auc_rsf_cause1[i, ] <- list_rsf[[i]]$AUC_cause1[1:11]
  auc_rsf_cause2[i, ] <- list_rsf[[i]]$AUC_cause2[1:11]
  
  calib_rsf2y_cause1[i] <- list_rsf[[i]]$MSE2_cause1
  calib_rsf2y_cause2[i] <- list_rsf[[i]]$MSE2_cause2
  calib_rsf5y_cause1[i] <- list_rsf[[i]]$MSE5_cause1
  calib_rsf5y_cause2[i] <- list_rsf[[i]]$MSE5_cause2
  calib_rsf10y_cause1[i] <- list_rsf[[i]]$MSE10_cause1
  calib_rsf10y_cause2[i] <- list_rsf[[i]]$MSE10_cause2 
  
}


# Boxplots for miscalibration (for the event of interest: disease progression)

# at 2 years

df1 <- data.frame(Calib_groups = calib_csc2y_cause1,
                  Model = "Cause-specific Cox",
                  Time = "2 years")


df2 <- data.frame(Calib_groups = calib_fg2y_cause1,
                  Model = "Fine-Gray",
                  Time = "2 years")

df3 <- data.frame(Calib_groups = calib_nnet2y_cause1,
                  Model = "PLANNCR original",
                  Time = "2 years")


df4 <- data.frame(Calib_groups = calib_keras2y_cause1,
                  Model = "PLANNCR extended",
                  Time = "2 years")

df5 <- data.frame(Calib_groups = calib_rsf2y_cause1,
                  Model = "RSFCR",
                  Time = "2 years")


#at 5 years

df6 <- data.frame(Calib_groups = calib_csc5y_cause1,
                  Model = "Cause-specific Cox",
                  Time = "5 years")


df7 <- data.frame(Calib_groups = calib_fg5y_cause1,
                  Model = "Fine-Gray",
                  Time = "5 years")

df8 <- data.frame(Calib_groups = calib_nnet5y_cause1,
                  Model = "PLANNCR original",
                  Time = "5 years")


df9 <- data.frame(Calib_groups = calib_keras5y_cause1,
                  Model = "PLANNCR extended",
                  Time = "5 years")

df10 <- data.frame(Calib_groups = calib_rsf5y_cause1,
                   Model = "RSFCR",
                   Time = "5 years")

#at 10 years

df11 <- data.frame(Calib_groups = calib_csc10y_cause1,
                   Model = "Cause-specific Cox",
                   Time = "10 years")


df12 <- data.frame(Calib_groups = calib_fg10y_cause1,
                   Model = "Fine-Gray",
                   Time = "10 years")

df13 <- data.frame(Calib_groups = calib_nnet10y_cause1,
                   Model = "PLANNCR original",
                   Time = "10 years")


df14 <- data.frame(Calib_groups = calib_keras10y_cause1,
                   Model = "PLANNCR extended",
                   Time = "10 years")

df15 <- data.frame(Calib_groups = calib_rsf10y_cause1,
                   Model = "RSFCR",
                   Time = "10 years")


dfa <- rbind(df1, df2, df3, df4, df5, df6, df7, df8,
             df9, df10, df11, df12, df13, df14, df15)

dfa$Model <- as.factor(dfa$Model)
dfa$Model <- factor(dfa$Model, levels = c("Cause-specific Cox", "Fine-Gray", "PLANNCR original",
                                          "PLANNCR extended", "RSFCR"))
dfa$Time <- as.factor(dfa$Time)
dfa$Time <- factor(dfa$Time, levels = c("2 years", "5 years", "10 years"))


calib_c1 <- ggplot(dfa, aes(x=Time, y=Calib_groups, fill = Model)) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_brewer(palette="PuBu") +  
  ylim(c(0, 0.03)) + 
  labs(title = " ",x = "Time since surgery",
       y = "Miscalibration (for disease progression)", fill = " ") +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top", 
        plot.title = element_text(hjust = 0.5)) # Remove legend by setting to "none" 


calib_c1

# Boxplots for miscalibration (for the competing event: death)

# at 2 years

df16 <- data.frame(Calib_groups = calib_csc2y_cause2,
                   Model = "Cause-specific Cox",
                   Time = "2 years")


df17 <- data.frame(Calib_groups = calib_fg2y_cause2,
                   Model = "Fine-Gray",
                   Time = "2 years")

df18 <- data.frame(Calib_groups = calib_nnet2y_cause2,
                   Model = "PLANNCR original",
                   Time = "2 years")


df19 <- data.frame(Calib_groups = calib_keras2y_cause2,
                   Model = "PLANNCR extended",
                   Time = "2 years")

df20 <- data.frame(Calib_groups = calib_rsf2y_cause2,
                   Model = "RSFCR",
                   Time = "2 years")


#at 5 years

df21 <- data.frame(Calib_groups = calib_csc5y_cause2,
                   Model = "Cause-specific Cox",
                   Time = "5 years")


df22 <- data.frame(Calib_groups = calib_fg5y_cause2,
                   Model = "Fine-Gray",
                   Time = "5 years")

df23 <- data.frame(Calib_groups = calib_nnet5y_cause2,
                   Model = "PLANNCR original",
                   Time = "5 years")


df24 <- data.frame(Calib_groups = calib_keras5y_cause2,
                   Model = "PLANNCR extended",
                   Time = "5 years")

df25 <- data.frame(Calib_groups = calib_rsf5y_cause2,
                   Model = "RSFCR",
                   Time = "5 years")

#at 10 years

df26 <- data.frame(Calib_groups = calib_csc10y_cause2,
                   Model = "Cause-specific Cox",
                   Time = "10 years")


df27 <- data.frame(Calib_groups = calib_fg10y_cause2,
                   Model = "Fine-Gray",
                   Time = "10 years")

df28 <- data.frame(Calib_groups = calib_nnet10y_cause2,
                   Model = "PLANNCR original",
                   Time = "10 years")


df29 <- data.frame(Calib_groups = calib_keras10y_cause2,
                   Model = "PLANNCR extended",
                   Time = "10 years")

df30 <- data.frame(Calib_groups = calib_rsf10y_cause2,
                   Model = "RSFCR",
                   Time = "10 years")


dfb <- rbind(df16, df17, df18, df19, df20, df21, df22, df23,
             df24, df25, df26, df27, df28, df29, df30)

dfb$Model <- as.factor(dfb$Model)
dfb$Model <- factor(dfb$Model, levels = c("Cause-specific Cox", "Fine-Gray", "PLANNCR original",
                                          "PLANNCR extended", "RSFCR"))
dfb$Time <- as.factor(dfb$Time)
dfb$Time <- factor(dfb$Time, levels = c("2 years", "5 years", "10 years"))



calib_c2 <- ggplot(dfb, aes(x=Time, y=Calib_groups, fill = Model)) + 
  geom_boxplot(notch = FALSE) +
  scale_fill_brewer(palette="OrRd") +  
  ylim(c(0, 0.03)) + 
  labs(title = " ",x = "Time since surgery",
       y = "Miscalibration (for death)", fill = " ") +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=11)) +
  theme(legend.position = "top", 
        plot.title = element_text(hjust = 0.5)) # Remove legend by setting to "none" 

calib_c2



## end of the file





