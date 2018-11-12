# R package conformal.glm 

## Conformal Prediction for Generalized Linear Regression Models

This package computes and compares prediction regions for the normal, Gamma, 
and inverse Gaussian families in the `glm` package.  There is 
functionality to construct the usual Wald type prediction region that one 
obtains from maximum likelihood estimation and the delta method, the 
parametric conformal prediction region, the nonparametric conformal 
prediction region, and prediction regions from conformalization of residuals. 


## Usage 

```r
library(devtools)
install_github(repo = "DEck13/conformal.glm", subdir="conformal.glm")
```

[simple example that illustrates functionality]
```r
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

fit = glm(y ~ x1, family = Gamma, data = data) 
cpred = conformal.glm(fit, nonparametric = TRUE, bins = 3, 
	cores = 6)
```

[least squares conformal prediction from the conformalInference package]
```r
library(conformalInference)
funs <- lm.funs(intercept = TRUE)
train.fun <- funs$train.fun
predict.fun <- funs$predict.fun
p1.tibs <- conformal.pred(x = x, y = y, x0 = x, 
  train.fun = train.fun, predict.fun = predict.fun, 
  alpha = alpha, grid.method = "linear",
  num.grid.pts = 999)
cresid = cbind(p1.tibs$lo, p1.tibs$up)
```

[plots of parametric, non parametric, and least squares precition regions and the plot of the prediction region obtained from the delta method]
```r
par(mfrow = c(2,2))
index <- sort(x, index.return = TRUE)$ix

## parametric conformal prediction region
paraCI <- cpred$paraconformal
mean(apply(paraCI, 1, diff))
plot(x, y, pch = 20)
lines(x[index], paraCI[, 1][index], type = "l", col = "red")
lines(x[index], paraCI[, 2][index], type = "l", col = "red")

## nonparametric conformal prediction region
nonparaCI <- cpred$nonparaconformal
mean(apply(nonparaCI, 1, diff))
plot(x, y, pch = 20)
lines(x[index], nonparaCI[, 1][index], type = "l", col = "red")
lines(x[index], nonparaCI[, 2][index], type = "l", col = "red")

## least squares conformal prediction region
mean(apply(cresid, 1, diff))
plot(x, y, pch = 20)
lines(x[index], cresid[, 1][index], type = "l", col = "red")
lines(x[index], cresid[, 2][index], type = "l", col = "red")

## delta method prediction region
p1 <- predict(fit, type = "response", se.fit = TRUE)
pred <- p1$fit
se <- p1$se.fit * sqrt(n)
deltaCI <- cbind(pred + se*qnorm(alpha/2), pred + se*qnorm(1-alpha/2))
mean(apply(deltaCI, 1, diff))
plot(x, y, pch = 20)
lines(x[index], deltaCI[, 1][index], type = "l", col = "red")
lines(x[index], deltaCI[, 2][index], type = "l", col = "red")
```

[approximate the highest density region]
```r
betaMLE <- coefficients(fit)
shapeMLE <- as.numeric(gamma.shape(fit)[1])
rateMLE <- cbind(1, x) %*% betaMLE * shapeMLE
foo <- cbind(rateMLE, shapeMLE)
bar <- seq(from = 1e-5, 
  to = qgamma(alpha - 1e-5, rate = foo[j,1], shape = foo[j,2]), 
  length = 100000)
p.lwr <- pgamma(bar, rate = foo[j,1], shape = foo[j,2])
p.upr <- p.lwr + 1 - alpha
y.lwr <- qgamma(p.lwr, rate = foo[j,1], shape = foo[j,2])
y.upr <- qgamma(p.upr, rate = foo[j,1], shape = foo[j,2])
min.index <- which.min(y.upr - y.lwr)
p.a <- p.lwr[min.index]
p.b <- p.upr[min.index]

minlengthCI <- do.call(rbind, lapply(1:nrow(x), FUN = function(j){
  a <- qgamma(p.a, rate = foo[j,1], shape = foo[j,2])
  b <- qgamma(p.b, rate = foo[j,1], shape = foo[j,2])
  c(a,b)
}))
```


[comparison of parametric conformal prediction region and highest density region]
```r
par(mfrow = c(1,2))

## parametric conformal prediction region
mean(apply(paraCI, 1, diff))
plot(x, y, pch = 20)
lines(x[index], paraCI[, 1][index], type = "l", col = "red")
lines(x[index], paraCI[, 2][index], type = "l", col = "red")

## minimum length prediction region under assumed model
mean(apply(minlengthCI, 1, diff))
plot(x, y, pch = 20)
lines(x[index], minlengthCI[, 1][index], type = "l", col = "red")
lines(x[index], minlengthCI[, 2][index], type = "l", col = "red")
```


To cite this package:
```r
citation("conformal.glm")
```


## Further details

For more details on the parametric conformal prediction region, see:

  Eck, D.J., Crawford, F.W., and Aronow, P.M. (2018+)
  Conformal prediction for exponential families and generalized linear models.
  Preprint available on request (email daniel.eck@yale.edu).

For more details on the prediciton region formed from conformalization of 
residuals, see:

  Lei, J., G'Sell, M., Rinaldo, A., Tibshirani, R., and Wasserman, L. (2016)
  Distribution-Free Predictive Inference for Regression. 
  https://arxiv.org/abs/1604.04173

For more details on the nonparametric conformal prediction region, see:

  Lei, J. and Wasserman, L. (2014)
  Distribution-Free Prediction Bands for Non-parametric Regression. 
  Journal of the Royal Statistical Society: Series B, 76(1), 71-96.

