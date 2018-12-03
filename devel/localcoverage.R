

library(MASS)
library(conformal.glm)
library(parallel)
library(conformalInference)

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
newdata <- matrix(seq(0.01, 0.99, by = 0.01), ncol = 1)
colnames(newdata) <- c("x1")

fit = glm(y ~ x1, family = Gamma, data = data) 


## estimate local coverage probabilities when the prediction 
## region is calculated at evenly spaced grid points.
system.time(cpred <- conformal.glm(fit, nonparametric = TRUE, bins = 3, 
	newdata = newdata, cores = 6))

funs <- lm.funs(intercept = TRUE)
train.fun <- funs$train.fun
predict.fun <- funs$predict.fun
system.time(p1.tibs <- conformal.pred(x = x, y = y, x0 = newdata, 
  train.fun = train.fun, predict.fun = predict.fun, 
  alpha = alpha, grid.method = "linear",
  num.grid.pts = 999))
cresid = cbind(p1.tibs$lo, p1.tibs$up)
cresid.region <- cbind(newdata, cresid)

k <- 1 # number of main effects
index.data <- find.index(x, wn = 1/3, k = k)
index.newdata <- find.index(newdata, wn = 1/3, k = k)
parametric.region <- cbind(newdata, cpred$paraconformal)
nonparametric.region <- cbind(newdata, cpred$nonparaconformal)

local.coverage(region = parametric.region, data = data, newdata = newdata, 
	bins = 3, k = 1,  at.data = FALSE)
local.coverage(region = nonparametric.region, data = data, newdata = newdata, 
	bins = 3, k = 1,  at.data = FALSE)
local.coverage(region = cresid.region, data = data, newdata = newdata, 
	bins = 3, k = 1,  at.data = FALSE)




## estimate local coverage probabilities at the observed data
system.time(cpred <- conformal.glm(fit, nonparametric = TRUE, 
	bins = 3, cores = 6))

index.data <- find.index(x, wn = 1/3, k = k)
index.newdata <- find.index(x, wn = 1/3, k = k)
parametric.region <- cbind(x, cpred$paraconformal)
nonparametric.region <- cbind(x, cpred$nonparaconformal)

system.time(p1.tibs <- conformal.pred(x = x, y = y, x0 = x, 
  train.fun = train.fun, predict.fun = predict.fun, 
  alpha = alpha, grid.method = "linear",
  num.grid.pts = 999))
cresid = cbind(p1.tibs$lo, p1.tibs$up)
cresid.region <- cbind(x, cresid)

local.coverage(region = parametric.region, data = data, newdata = x, 
	bins = 3, k = 1,  at.data = TRUE)
local.coverage(region = nonparametric.region, data = data, newdata = x, 
	bins = 3, k = 1,  at.data = TRUE)
local.coverage(region = cresid.region, data = data, newdata = x, 
	bins = 3, k = 1,  at.data = TRUE)


