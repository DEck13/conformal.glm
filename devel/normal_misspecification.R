

library(MASS)
library(conformal.glm)
library(parallel)

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

fit = glm(y ~ x1, family = gaussian, data = data) 
system.time(cpred <- conformal.glm(fit, nonparametric = TRUE, bins = 3, 
	family = "gaussian", newdata = newdata, cores = 6))



library(conformalInference)
funs <- lm.funs(intercept = TRUE)
train.fun <- funs$train.fun
predict.fun <- funs$predict.fun
system.time(p1.tibs <- conformal.pred(x = x, y = y, x0 = newdata, 
  train.fun = train.fun, predict.fun = predict.fun, 
  alpha = alpha, grid.method = "linear",
  num.grid.pts = 999))
cresid = cbind(p1.tibs$lo, p1.tibs$up)


par(mfrow = c(2,2))

## parametric conformal prediction region
paraCI <- cpred$paraconformal
mean(apply(paraCI, 1, diff))
plot(x, y, pch = 20)
lines(newdata, paraCI[, 1], type = "l", col = "red")
lines(newdata, paraCI[, 2], type = "l", col = "red")

## nonparametric conformal prediction region
nonparaCI <- cpred$nonparaconformal
mean(apply(nonparaCI, 1, diff))
plot(x, y, pch = 20)
lines(newdata, nonparaCI[, 1], type = "l", col = "red")
lines(newdata, nonparaCI[, 2], type = "l", col = "red")

## least squares conformal prediction region
mean(apply(cresid, 1, diff))
plot(x, y, pch = 20)
lines(newdata, cresid[, 1], type = "l", col = "red")
lines(newdata, cresid[, 2], type = "l", col = "red")

