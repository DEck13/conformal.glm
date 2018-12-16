

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

## parametric and nonparametric conformal prediction region
system.time(cpred <- conformal.glm(fit, nonparametric = TRUE, bins = 3, 
	newdata = newdata, cores = 6))

## least squares conformal prediction region
library(conformalInference)
funs <- lm.funs(intercept = TRUE)
train.fun <- funs$train.fun
predict.fun <- funs$predict.fun
system.time(p1.tibs <- conformal.pred(x = x, y = y, x0 = newdata, 
  train.fun = train.fun, predict.fun = predict.fun, 
  alpha = alpha, grid.method = "linear",
  num.grid.pts = 999))
cresid = cbind(p1.tibs$lo, p1.tibs$up)

## highest density region
library(HDInterval)
betaMLE <- coefficients(fit)
shapeMLE <- as.numeric(gamma.shape(fit)[1])
rateMLE <- cbind(1, newdata) %*% betaMLE * shapeMLE
minlength <- do.call(rbind, 
  lapply(1:nrow(newdata), function(j){ 
    hdi(qgamma, 0.90, shape = shapeMLE, rate = rateMLE[j, 1])
  }))



#########################################
## area and coverage

## parametric conformal prediction region
paraCI <- cpred$paraconformal
# estimated area
mean(apply(paraCI, 1, diff))
# local coverage
p <- length(beta) - 1
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


#########################################
## make plot
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




file <- paste("readme", collapse = "")
file <- paste(file, "RData", sep = ".", collapse = "")
save.image(file = file, ascii = TRUE)

