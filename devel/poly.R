


library(conformal.glm)
library(statmod)
library(MASS)
library(parallel)

set.seed(13)
n <- 200
shape <- 2
beta <- c(1, 1, 1, 1)
x1 <- runif(n)
x2 <- runif(n)
rate <- cbind(1, x1, x2, x1^2) %*% beta * shape
y <- rgamma(n = n, shape = shape, rate = rate)
data <- data.frame(y = y, x1 = x1, x2 = x2)
fit = glm(y ~ x1 + x2 + I(x1^2), family = "Gamma", data = data) 

gout <- conformal.glm(fit, bins = 2, cores = 6, 
  nonparametric = TRUE)
gout
