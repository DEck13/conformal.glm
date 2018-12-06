

library(parallel)
library(MASS)
library(statmod)
library(conformal.glm)


alpha <- 0.10
n <- 200
shape <- 2
beta <- c(1/4, 2)

set.seed(13)
x <- matrix(runif(n), ncol = 1)
rate <- cbind(1, x) %*% beta * shape
y <- rgamma(n = n, shape = shape, rate = rate)
data <- data.frame(y = y, x = x)
colnames(data)[2] <- c("x1")

fit = glm(y ~ x1, family = Gamma, data = data) 
cpred <- conformal.glm(fit, nonparametric = TRUE, 
	bins = 3, cores = 6)
paraCI <- cpred$paraconformal



## Monte Carlo coverage 
montecarlo <- function(n, bins, beta = beta, shape = shape){
	p <- length(beta) - 1
  x <- matrix(runif(n), ncol = p)
  rate <- cbind(1, x) %*% beta * shape
  y <- rgamma(n = n, shape = shape, rate = rate)
  data <- data.frame(y = y, x = x)
  colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")

  fit = glm(y ~ ., family = Gamma, data = data) 
  cpred <- conformal.glm(fit, nonparametric = FALSE, 
	  bins = bins, cores = 6)
  paraCI <- cpred$paraconformal
  nonparaCI <- cpred$nonparaconformal

  formula <- fit$formula
  newdata <- data
  respname <- all.vars(formula)[1]
  newdata <- newdata[, !(colnames(data) %in% respname)]
  newdata <- as.matrix(newdata)

  ## local coverage for parametric conformal prediction region
  local.para <- local.coverage(region = cbind(newdata, paraCI), 
  	data = data, newdata = newdata, k = p, bins = bins, at.data = "TRUE")
  c(local.para, mean(apply(paraCI, 1, diff)))

}






## n = 200, bins = 3
B <- 25
n <- 200
bins <- 3
system.time(mat200.3 <- matrix(unlist(lapply(1:B, 
	FUN = function(j) montecarlo(n = n, bins = bins, beta = beta, shape = shape))), 
  byrow = T, nrow = B))
apply(mat200.3, 2, mean)
apply(mat200.3, 2, sd)

## n = 200, bins = 4
B <- 25
n <- 200
bins <- 4
system.time(mat200.4 <- matrix(unlist(lapply(1:B, 
	FUN = function(j) montecarlo(n = n, bins = bins, beta = beta, shape = shape))), 
  byrow = T, nrow = B))
apply(mat200.4, 2, mean)
apply(mat200.4, 2, sd)



## n = 500, bins = 4
B <- 25
n <- 500
bins <- 4
system.time(mat500.4 <- matrix(unlist(lapply(1:B, 
	FUN = function(j) montecarlo(n = n, bins = bins, beta = beta, shape = shape))), 
  byrow = T, nrow = B))
apply(mat500.4, 2, mean)
apply(mat500.4, 2, sd)

## n = 500, bins = 6
B <- 25
n <- 500
bins <- 6
system.time(mat500.6 <- matrix(unlist(lapply(1:B, 
	FUN = function(j) montecarlo(n = n, bins = bins, beta = beta, shape = shape))), 
  byrow = T, nrow = B))
apply(mat500.6, 2, mean)
apply(mat500.6, 2, sd)



## n = 1000, bins = 6
B <- 25
n <- 1000
bins <- 6
system.time(mat1000.6 <- matrix(unlist(lapply(1:B, 
	FUN = function(j) montecarlo(n = n, bins = bins, beta = beta, shape = shape))), 
  byrow = T, nrow = B))
apply(mat1000.6, 2, mean)
apply(mat1000.6, 2, sd)

## n = 1000, bins = 8
B <- 25
n <- 1000
bins <- 8
system.time(mat1000.8 <- matrix(unlist(lapply(1:B, 
	FUN = function(j) montecarlo(n = n, bins = bins, beta = beta, shape = shape))), 
  byrow = T, nrow = B))
apply(mat1000.8, 2, mean)
apply(mat1000.8, 2, sd)



## n = 2000, bins = 8
B <- 25
n <- 2000
bins <- 8
system.time(mat2000.8 <- matrix(unlist(lapply(1:B, 
	FUN = function(j) montecarlo(n = n, bins = bins, beta = beta, shape = shape))), 
  byrow = T, nrow = B))
apply(mat2000.8, 2, mean)
apply(mat2000.8, 2, sd)

## n = 2000, bins = 10
B <- 25
n <- 2000
bins <- 10
system.time(mat2000.10 <- matrix(unlist(lapply(1:B, 
	FUN = function(j) montecarlo(n = n, bins = bins, beta = beta, shape = shape))), 
  byrow = T, nrow = B))
apply(mat2000.10, 2, mean)
apply(mat2000.10, 2, sd)

## n = 2000, bins = 12
B <- 25
n <- 2000
bins <- 12
system.time(mat2000.12 <- matrix(unlist(lapply(1:B, 
	FUN = function(j) montecarlo(n = n, bins = bins, beta = beta, shape = shape))), 
  byrow = T, nrow = B))
apply(mat2000.12, 2, mean)
apply(mat2000.12, 2, sd)




file <- paste("Gamma-parametric-simulation", collapse = "")
file <- paste(file, "RData", sep = ".", collapse = "")
save.image(file = file, ascii = TRUE)


