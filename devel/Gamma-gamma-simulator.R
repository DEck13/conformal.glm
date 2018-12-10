

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
shape <- 2
beta <- c(1/8, 2)



## n = 100, bins = 2
B <- 50
n <- 100
bins <- 2
system.time(gamgam100.2 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    print(j)
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "Gamma", 
      link = "inverse", shape = shape, sd = NULL, 
      confamily = "Gamma",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gamgam100.2[, 1])))
B <- B - NAs
cbind(apply(gamgam100.2, 2, function(j) mean(j, na.rm = TRUE) ), 
  apply(gamgam100.2, 2, function(j) sd(j, na.rm = TRUE)) / sqrt(B))


## n = 200, bins = 2
B <- 50
n <- 200
bins <- 2
system.time(gamgam200.2 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "Gamma", 
      link = "inverse", shape = shape, sd = NULL, 
      confamily = "Gamma",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gamgam200.2[, 1])))
B <- B - NAs
cbind(apply(gamgam200.2, 2, mean), 
  apply(gamgam200.2, 2, sd) / sqrt(B))


## n = 200, bins = 3
B <- 50
n <- 200
bins <- 3
system.time(gamgam200.3 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "Gamma", 
      link = "inverse", shape = shape, sd = NULL, 
      confamily = "Gamma",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gamgam200.3[, 1])))
B <- B - NAs
cbind(apply(gamgam200.3, 2, mean), 
  apply(gamgam200.3, 2, sd) / sqrt(B))


## n = 500, bins = 2
B <- 50
n <- 500
bins <- 2
system.time(gamgam500.2 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "Gamma", 
      link = "inverse", shape = shape, sd = NULL, 
      confamily = "Gamma",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gamgam500.2[, 1])))
B <- B - NAs
cbind(apply(gamgam500.2, 2, mean), 
  apply(gamgam500.2, 2, sd) / sqrt(B))


## n = 500, bins = 3
B <- 50
n <- 500
bins <- 3
system.time(gamgam500.3 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "Gamma", 
      link = "inverse", shape = shape, sd = NULL, 
      confamily = "Gamma",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gamgam500.3[, 1])))
B <- B - NAs
cbind(apply(gamgam500.3, 2, mean), 
  apply(gamgam500.3, 2, sd) / sqrt(B))


## n = 500, bins = 4
B <- 50
n <- 500
bins <- 4
system.time(gamgam500.4 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "Gamma", 
      link = "inverse", shape = shape, sd = NULL, 
      confamily = "Gamma",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gamgam500.4[, 1])))
B <- B - NAs
cbind(apply(gamgam500.4, 2, mean), 
  apply(gamgam500.4, 2, sd) / sqrt(B))


## n = 500, bins = 5
B <- 50
n <- 500
bins <- 5
system.time(gamgam500.5 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    unlist(simulator(n = n, alpha = alpha, beta = beta, 
      bins = bins, family = "Gamma", 
      link = "inverse", shape = shape, sd = NULL, 
      confamily = "Gamma",
      parametric = TRUE, nonparametric = TRUE, 
      LS = TRUE, HD = TRUE, cores = 6))
  })))
NAs <- length(which(is.na(gamgam500.5[, 1])))
B <- B - NAs
cbind(apply(gamgam500.5, 2, mean), 
  apply(gamgam500.5, 2, sd) / sqrt(B))







file <- paste("Gamma-gamma-output", collapse = "")
file <- paste(file, "RData", sep = ".", collapse = "")
save.image(file = file, ascii = TRUE)


