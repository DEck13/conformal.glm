

library(parallel)
library(MASS)
library(statmod)
library(conformal.glm)
library(conformalInference)
library(HDInterval)
source("simulator.R")
set.seed(13)



## inputs
alpha <- 0.10
beta <- c(1, 4)
sd <- 2


## n = 100, bins = 2
B <- 50
n <- 100
bins <- 2
system.time(gauss100.2 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    print(j)
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "gaussian", 
      link = "inverse", shape = NULL, sd = sd, 
      confamily = "gaussian",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gauss100.2[, 1])))
B <- B - NAs
cbind(apply(gauss100.2, 2, function(j) mean(j, na.rm = TRUE) ), 
  apply(gauss100.2, 2, function(j) sd(j, na.rm = TRUE)) / sqrt(B))


## n = 200, bins = 2
B <- 50
n <- 200
bins <- 2
system.time(gauss200.2 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "gaussian", 
      link = "inverse", shape = NULL, sd = sd, 
      confamily = "gaussian",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gauss200.2[, 1])))
B <- B - NAs
cbind(apply(gauss200.2, 2, mean), 
  apply(gauss200.2, 2, sd) / sqrt(B))


## n = 250, bins = 2
B <- 50
n <- 250
bins <- 2
system.time(gauss250.2 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "gaussian", 
      link = "inverse", shape = NULL, sd = sd, 
      confamily = "gaussian",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gauss250.2[, 1])))
B <- B - NAs
cbind(apply(gauss250.2, 2, mean), 
  apply(gauss250.2, 2, sd) / sqrt(B))


## n = 300, bins = 2
B <- 50
n <- 300
bins <- 2
system.time(gauss300.2 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "gaussian", 
      link = "inverse", shape = NULL, sd = sd, 
      confamily = "gaussian",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gauss300.2[, 1])))
B <- B - NAs
cbind(apply(gauss300.2, 2, mean), 
  apply(gauss300.2, 2, sd) / sqrt(B))


## n = 400, bins = 2
B <- 50
n <- 400
bins <- 2
system.time(gauss400.2 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "gaussian", 
      link = "inverse", shape = NULL, sd = sd, 
      confamily = "gaussian",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gauss400.2[, 1])))
B <- B - NAs
cbind(apply(gauss400.2, 2, mean), 
  apply(gauss400.2, 2, sd) / sqrt(B))


## n = 400, bins = 3
B <- 50
n <- 400
bins <- 3
system.time(gauss400.3 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "gaussian", 
      link = "inverse", shape = NULL, sd = sd, 
      confamily = "gaussian",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gauss400.3[, 1])))
B <- B - NAs
cbind(apply(gauss400.3, 2, mean), 
  apply(gauss400.3, 2, sd) / sqrt(B))


## n = 500, bins = 2
B <- 50
n <- 500
bins <- 2
system.time(gauss500.2 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "gaussian", 
      link = "inverse", shape = NULL, sd = sd, 
      confamily = "gaussian",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gauss500.2[, 1])))
B <- B - NAs
cbind(apply(gauss500.2, 2, mean), 
  apply(gauss500.2, 2, sd) / sqrt(B))


## n = 500, bins = 3
B <- 50
n <- 500
bins <- 3
system.time(gauss500.3 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "gaussian", 
      link = "inverse", shape = NULL, sd = sd, 
      confamily = "gaussian",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gauss500.3[, 1])))
B <- B - NAs
cbind(apply(gauss500.3, 2, mean), 
  apply(gauss500.3, 2, sd) / sqrt(B))






file <- paste("gaussian-output", collapse = "")
file <- paste(file, "RData", sep = ".", collapse = "")
save.image(file = file, ascii = TRUE)

