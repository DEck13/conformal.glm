

## Test for speed: (n = 2000, bins = 8), 
# first run: 53.610 seconds with old code 
# conservative region: trimmed down to ~18 seconds!
# double back line search: 44.578 seconds (~45 seconds)

library(MASS)
library(conformal.glm)
library(parallel)

alpha <- 0.10
n <- 2000
shape <- 2
beta <- c(1/4, 2)

set.seed(13)
x <- matrix(runif(n), ncol = 1)
rate <- cbind(1, x) %*% beta * shape
y <- rgamma(n = n, shape = shape, rate = rate)
data <- data.frame(y = y, x = x)
colnames(data)[2] <- c("x1")
newdata <- matrix(seq(0.01, 0.99, by = 0.01), ncol = 1)
colnames(newdata) <- c("x1")

fit = glm(y ~ x1, family = Gamma, data = data) 
system.time(cpred <- conformal.glm(fit, nonparametric = FALSE, bins = 8, 
	newdata = newdata, cores = 6))

system.time(cpred <- conformal.glm(fit, parametric = FALSE, 
	nonparametric = TRUE, bins = 8, 
	newdata = newdata, cores = 6))

