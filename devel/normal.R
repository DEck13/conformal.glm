

library(parallel)
library(conformal.glm)


# ---- generate X data -------
set.seed(13)
alpha <- 0.10
n <- 5e2
beta <- c(1, 10)
sd <- 1/2
x <- matrix(runif(n), ncol = 1)
p <- ncol(x)
colnames(x) <- "x1"

# ---- generate regression data -------
y <- rnorm(n = n, mean = cbind(1 , x) %*% beta, sd = sd)
data <- data.frame(y = y, x = x)
fit <- glm(y ~ x1, data = data, family = "gaussian")


# ---- get new data for predictions -------
newdata <- cbind(seq(from = 0.01, to = 0.99, by = 0.01))
colnames(newdata) <- "x1"

# ---- parametric conformal prediction region -------
system.time(out <- conformal.glm(fit, 
	newdata = newdata, bins = 5))
paraCI <- out$paraconformal

# ---- nonparametric conformal prediction region -------
system.time(out2 <- conformal.glm(fit, newdata = newdata, 
	bins = 5, parametric = FALSE, nonparametric = TRUE))
nonparaCI <- out2$nonparaconformal

# ---- LS conformal prediction region -------
library(conformalInference)
funs <- lm.funs(intercept = TRUE)
train.fun <- funs$train.fun
predict.fun <- funs$predict.fun
system.time(p1.tibs <- conformal.pred(x = x, y = y, 
	x0 = newdata, train.fun = train.fun, 
	predict.fun = predict.fun, 
  alpha = alpha, grid.method = "linear",
  num.grid.pts = 999))
cresid = cbind(p1.tibs$lo, p1.tibs$up)



## parametric conformal prediction region
par(mfrow = c(2,2))

mean(apply(paraCI, 1, diff))
plot(x, y, pch = 20, main = "parametric conformal")
lines(newdata, paraCI[, 1], type = "l", col = "red")
lines(newdata, paraCI[, 2], type = "l", col = "red")

## nonparametric conformal prediction region
mean(apply(nonparaCI, 1, diff))
plot(x, y, pch = 20, main = "nonparametric conformal")
lines(newdata, nonparaCI[, 1], type = "l", col = "red")
lines(newdata, nonparaCI[, 2], type = "l", col = "red")

## nonparametric conformal prediction region
mean(apply(cresid, 1, diff))
plot(x, y, pch = 20, main = "LS conformal")
lines(newdata, cresid[, 1], type = "l", col = "red")
lines(newdata, cresid[, 2], type = "l", col = "red")
