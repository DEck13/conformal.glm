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

## Illustrative Example 

We provide a gamma regression example with perfect model specification to 
illustrate the performance of conformal predictions when the model is known 
and the model does not have additive symmetric errors.  We also compare 
conformal prediction regions to the oracle highest density region under 
the correct model. This example is included in the corresponding paper:  

  Eck, D.J., Crawford, F.W., and Aronow, P.M. (2018+)
  Conformal prediction for exponential families and generalized linear models.
  Preprint available on request (email daniel.eck@yale.edu)..

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
```


# parametric and nonparametric conformal prediction regions
```r
system.time(cpred <- conformal.glm(fit, nonparametric = TRUE, bins = 3, 
  newdata = newdata, cores = 6))
```


# least squares conformal prediction from the conformalInference package 
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


# estimate the highest density region
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


## Plot of all prediction regions
```r
par(mfrow = c(2,2), oma = c(4,4,0,0), mar = c(1,1,1,1))

# parametric conformal prediction region
plot(x, y, pch = 20, xaxt = 'n', ann = FALSE)
lines(newdata, paraCI[, 1], type = "l", col = "red")
lines(newdata, paraCI[, 2], type = "l", col = "red")

# nonparametric conformal prediction region
plot(x, y, pch = 20, xaxt = 'n', yaxt = 'n', ann = FALSE)
lines(newdata, nonparaCI[, 1], type = "l", col = "red")
lines(newdata, nonparaCI[, 2], type = "l", col = "red")

# least squares conformal prediction region
plot(x, y, pch = 20, ann = FALSE)
lines(newdata, cresid[, 1], type = "l", col = "red")
lines(newdata, cresid[, 2], type = "l", col = "red")

# highest density region
plot(x, y, pch = 20, yaxt = 'n', ann = FALSE)
lines(newdata, minlength[, 1], type = "l", col = "red")
lines(newdata, minlength[, 2], type = "l", col = "red")

# axis labels
mtext("x", side = 1, line = 2.5, outer = TRUE, cex = 2)
mtext("y", side = 2, line = 2.5, outer = TRUE, cex = 2)
```

![Depiction of prediction regions](https://github.com/DEck13/conformal.glm/tree/master/gammasimexample.pdf)


## Coverage properties and estimated area of all prediction regions
```r
## parametric conformal prediction region
paraCI <- cpred$paraconformal
# estimated area
mean(apply(paraCI, 1, diff))
# local coverage
local.coverage(region = paraCI, 
      data = data, newdata = newdata, k = p, bins = 3, 
      at.data = "FALSE")
# marginal coverage
local.coverage(region = paraCI, 
      data = data, newdata = newdata, k = p, bins = 1, 
      at.data = "FALSE")

## nonparametric conformal prediction region
nonparaCI <- cpred$nonparaconformal
# estimated area
mean(apply(nonparaCI, 1, diff))
# local coverage
local.coverage(region = nonparaCI, 
      data = data, newdata = newdata, k = p, bins = 3, 
      at.data = "FALSE")
# marginal coverage
local.coverage(region = nonparaCI, 
      data = data, newdata = newdata, k = p, bins = 1, 
      at.data = "FALSE")

## least squares conformal prediction region
# estimated area
mean(apply(cresid, 1, diff))
# local coverage
local.coverage(region = cresid, 
      data = data, newdata = newdata, k = p, bins = 3, 
      at.data = "FALSE")
# marginal coverage
local.coverage(region = cresid, 
      data = data, newdata = newdata, k = p, bins = 1, 
      at.data = "FALSE")

## highest density region
# estimated area
mean(apply(minlength, 1, diff))
# local coverage
local.coverage(region = minlength, 
      data = data, newdata = newdata, k = p, bins = 3, 
      at.data = "FALSE")
# marginal coverage
local.coverage(region = minlength, 
      data = data, newdata = newdata, k = p, bins = 1, 
      at.data = "FALSE")
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

