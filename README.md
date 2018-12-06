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
newdata <- matrix(seq(0.01, 0.99, by = 0.01), ncol = 1)
colnames(newdata) <- c("x1")

fit = glm(y ~ x1, family = Gamma, data = data) 
system.time(cpred <- conformal.glm(fit, nonparametric = TRUE, bins = 3, 
	newdata = newdata, cores = 6))
```

[least squares conformal prediction from the conformalInference package]
```r
library(conformalInference)
funs <- lm.funs(intercept = TRUE)
train.fun <- funs$train.fun
predict.fun <- funs$predict.fun
system.time(p1.tibs <- conformal.pred(x = x, y = y, x0 = newdata, 
  train.fun = train.fun, predict.fun = predict.fun, 
  alpha = alpha, grid.method = "linear",
  num.grid.pts = 999))
cresid = cbind(p1.tibs$lo, p1.tibs$up)
```

[plots of parametric, non parametric, and least squares conformal precition regions]
```r
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
```

[estimate the highest density region]
```r
library(HDInterval)
betaMLE <- coefficients(fit)
shapeMLE <- as.numeric(gamma.shape(fit)[1])
rateMLE <- cbind(1, newdata) %*% betaMLE * shapeMLE

minlength <- do.call(rbind, 
  lapply(1:nrow(newdata), function(j){ 
    hdi(qgamma, 0.90, shape = shapeMLE, rate = rateMLE[j, 1])
  }))
```


[plot the highest density region]
```r
## minimum length prediction region under assumed model
mean(apply(minlength, 1, diff))
plot(x, y, pch = 20)
lines(newdata, minlength[, 1], type = "l", col = "red")
lines(newdata, minlength[, 2], type = "l", col = "red")
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

