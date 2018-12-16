

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
beta <- c(1, 1/2)
sd <- 2


#par(mfrow = c(5,5))
#for(j in 1:length(shape)){
#  hist(rgamma(n=n, shape=shape[j], rate=rate))
#}

misgamgauss <- function(shape = shape, bins = bins){
  unlist(simulator(n = n, alpha = alpha, beta = beta, 
    bins = bins, family = "Gamma", 
    link = "inverse", shape = shape,
    confamily = "gaussian",
    parametric = TRUE, nonparametric = TRUE, 
    LS = TRUE, HD = TRUE, cores = 6))
}


## n = 500, bins = 2, shape = 2
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.2 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 2, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.2[, 1])))
B <- B - NAs
out500.2.2 <- cbind(apply(gamgauss500.2.2, 2, mean), 
  apply(gamgauss500.2.2, 2, sd) / sqrt(B))
rownames(out500.2.2) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.2



## n = 500, bins = 2, shape = 4
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.4 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 4, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.4[, 1])))
B <- B - NAs
out500.2.4 <- cbind(apply(gamgauss500.2.4, 2, mean), 
  apply(gamgauss500.2.4, 2, sd) / sqrt(B))
rownames(out500.2.4) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.4



## n = 500, bins = 2, shape = 6
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.6 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 6, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.6[, 1])))
B <- B - NAs
out500.2.6 <- cbind(apply(gamgauss500.2.6, 2, mean), 
  apply(gamgauss500.2.6, 2, sd) / sqrt(B))
rownames(out500.2.6) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.6



## n = 500, bins = 2, shape = 8
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.8 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 8, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.8[, 1])))
B <- B - NAs
out500.2.8 <- cbind(apply(gamgauss500.2.8, 2, mean), 
  apply(gamgauss500.2.8, 2, sd) / sqrt(B))
rownames(out500.2.8) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.8



## n = 500, bins = 2, shape = 10
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.10 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 10, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.10[, 1])))
B <- B - NAs
out500.2.10 <- cbind(apply(gamgauss500.2.10, 2, mean), 
  apply(gamgauss500.2.10, 2, sd) / sqrt(B))
rownames(out500.2.10) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.10



## n = 500, bins = 2, shape = 12
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.12 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 12, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.12[, 1])))
B <- B - NAs
out500.2.12 <- cbind(apply(gamgauss500.2.12, 2, mean), 
  apply(gamgauss500.2.12, 2, sd) / sqrt(B))
rownames(out500.2.12) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.12



## n = 500, bins = 2, shape = 14
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.14 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 14, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.14[, 1])))
B <- B - NAs
out500.2.14 <- cbind(apply(gamgauss500.2.14, 2, mean), 
  apply(gamgauss500.2.14, 2, sd) / sqrt(B))
rownames(out500.2.14) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.14



## n = 500, bins = 2, shape = 16
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.16 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 16, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.16[, 1])))
B <- B - NAs
out500.2.16 <- cbind(apply(gamgauss500.2.16, 2, mean), 
  apply(gamgauss500.2.16, 2, sd) / sqrt(B))
rownames(out500.2.16) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.16



## n = 500, bins = 2, shape = 18
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.18 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 18, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.18[, 1])))
B <- B - NAs
out500.2.18 <- cbind(apply(gamgauss500.2.18, 2, mean), 
  apply(gamgauss500.2.18, 2, sd) / sqrt(B))
rownames(out500.2.18) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.18



## n = 500, bins = 2, shape = 20
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.20 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 20, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.20[, 1])))
B <- B - NAs
out500.2.20 <- cbind(apply(gamgauss500.2.20, 2, mean), 
  apply(gamgauss500.2.20, 2, sd) / sqrt(B))
rownames(out500.2.20) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.20



file <- paste("misspecification-gamma-gaussian-output", collapse = "")
file <- paste(file, "RData", sep = ".", collapse = "")
save.image(file = file, ascii = TRUE)

