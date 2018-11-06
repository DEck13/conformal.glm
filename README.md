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
set.seed(13)
n <- 200
shape <- 1
beta <- c(1, 1/4)
x <- matrix(runif(n), ncol = 1)
rate <- cbind(1, x) %*% beta * shape
y <- rgamma(n = n, shape = shape, rate = rate)
data <- data.frame(y = y, x = x)
colnames(data)[2] <- c("x1")

fit = glm(y ~ x1, family = "Gamma", data = data) 
cpred = conformal.glm(fit, nonparametric = TRUE, bins = 5, 
	cores = 6)
```
[then show a picture or something]

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

