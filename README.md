# R package conformal.glm 

## Conformal Prediction for Generalized Linear Regression Models

This package computes and compares prediction regions for the normal, Gamma, 
and inverse Gaussian families in the `glm` package.  There is functionality to 
construct the parametric conformal prediction region and the nonparametric 
conformal prediction region. 


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
  Preprint available on request (email daniel.eck@yale.edu).

```r
library(MASS)
library(conformal.glm)
library(parallel)

alpha <- 0.10
n <- 500
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


Compute the parametric and nonparametric conformal prediction regions.
```r
system.time(cpred <- conformal.glm(fit, nonparametric = TRUE, bins = 5, 
  newdata = newdata, cores = 6))
paraCI <- cpred$paraconformal
nonparaCI <- cpred$nonparaconformal
```


Compute the least squares conformal prediction region from the conformalInference package.
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


Compute the highest density region.
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

The four prediction regions for this data are depicted below.  The top row 
depicts the parametric conformal prediction region (left panel) and the 
nonparametric conformal prediction region (right panel).  The bin width was 
specified as 1/3 for these conformal prediction regions.  The bottom row 
depicts the least squares conformal prediction region (left panel) and the 
highest density region (right panel).  We see that the parametric conformal 
prediction region is a close discretization of the highest density region, the 
nonparametric conformal prediction region is quite jagged and unnatural, and 
the least squares conformal prediction region exhibits undercoverage for small 
x, exhibits overcoverage for large x, and includes negative response values 
of large magnitude. 

```r
par(mfrow = c(2,2), oma = c(4,4,0,0), mar = c(1,1,1,1))

# parametric conformal prediction region
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,max(y)))
points(x, y, pch = 19, col = "gray")
lines(newdata, paraCI[, 1], type = "l", col = "red")
lines(newdata, paraCI[, 2], type = "l", col = "red")
axis(2)

# nonparametric conformal prediction region
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,max(y)))
points(x, y, pch = 19, col = "gray")
lines(newdata, nonparaCI[, 1], type = "l", col = "red")
lines(newdata, nonparaCI[, 2], type = "l", col = "red")

# least squares conformal prediction region
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,max(y)))
points(x, y, pch = 19, col = "gray")
lines(newdata, cresid[, 1], type = "l", col = "red")
lines(newdata, cresid[, 2], type = "l", col = "red")
axis(1); axis(2)

# highest density region
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,max(y)))
points(x, y, pch = 19, col = "gray")
lines(newdata, minlength[, 1], type = "l", col = "red")
lines(newdata, minlength[, 2], type = "l", col = "red")
axis(1)

# axis labels
mtext("x", side = 1, line = 2.5, outer = TRUE, cex = 2)
mtext("y", side = 2, line = 2.5, outer = TRUE, cex = 2)
```


![Plot of prediction regions](gammasimexample.png)


## Coverage properties and estimated area of all prediction regions

All of the presented prediction regions exhibit finite-sample marginal 
validity.  However, the least squares conformal prediction region does not 
exhibit finite-sample local validity with bins of (0,1/3], (1/3, 2/3], 
(2/3, 1] used to assess local validity.  The highest density prediction 
region is the smallest in size with an estimated area of 2.28.  The parametric 
conformal prediction region is close in size with an estimated area of 2.35.  
The nonparametric conformal prediction region has an estimated area of 2.62 
and the least square conformal prediction region has an estimated area of 
2.94.  Under correct model specification, the parametric conformal prediction 
region is similar in performance to that of the highest density prediction 
region.


```r
## parametric conformal prediction region
paraCI <- cpred$paraconformal
# estimated area
mean(apply(paraCI, 1, diff))
# local coverage
p <- length(beta) - 1
local.coverage(region = paraCI, 
      data = data, newdata = newdata, k = p, bins = 5, 
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
      data = data, newdata = newdata, k = p, bins = 5, 
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
      data = data, newdata = newdata, k = p, bins = 5, 
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
      data = data, newdata = newdata, k = p, bins = 5, 
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

