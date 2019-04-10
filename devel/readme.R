

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
data.readme <- data.frame(y = y, x = x)
colnames(data.readme)[2] <- c("x1")

fit.readme = glm(y ~ x1, family = Gamma, data = data.readme)
summary(fit.readme)


## parametric and nonparametric conformal prediction region
bins <- 3
system.time(cpred <- conformal.glm(fit.readme, 
  nonparametric = TRUE, bins = bins, cores = 6))
paraCI <- cpred$paraconformal
nonparaCI <- cpred$nonparaconformal


## least squares conformal prediction region
library(conformalInference)
funs <- lm.funs(intercept = TRUE)
train.fun <- funs$train.fun
predict.fun <- funs$predict.fun

cubic.model <- lm(y ~ x + I(x^2) + I(x^3))
abs.resid <- abs(cubic.model$resid)
smooth.call <- smooth.spline(x, abs.resid, 
  nknots = 10)
lambda <- smooth.call$lambda
df <- smooth.call$df
mad.train.fun <- function(x, y, out = NULL){
  smooth.spline(x[, 1], y, lambda = lambda, 
    df = df, nknots = 10)
}
mad.predict.fun <- function(out, newx){
  predict(out, newx[, 1])$y
}
system.time(p1.tibs <- conformal.pred(x = cbind(x,x^2,x^3), 
  y = y, x0 = cbind(x,x^2,x^3), 
  train.fun = train.fun, predict.fun = predict.fun, 
  mad.train.fun = mad.train.fun,
  mad.predict.fun = mad.predict.fun,
  alpha = alpha))
LSLW = cbind(p1.tibs$lo, p1.tibs$up)


## highest density region
library(HDInterval)
betaMLE <- coefficients(fit.readme)
shapeMLE <- as.numeric(gamma.shape(fit.readme)[1])
rateMLE <- cbind(1, x) %*% betaMLE * shapeMLE
HDCI <- do.call(rbind, 
  lapply(1:nrow(x), function(j){ 
    hdi(qgamma, 0.90, shape = shapeMLE, rate = rateMLE[j, 1])
  }))



#########################################
## area and coverage

## parametric conformal prediction region
# estimated area
mean(apply(paraCI, 1, diff))
# local coverage
p <- length(beta) - 1
local.coverage(region = paraCI, data = data.readme, d = p, 
  bins = bins, at.data = "TRUE")
# marginal coverage
local.coverage(region = paraCI, data = data.readme, d = p, 
  bins = 1, at.data = "TRUE")

## nonparametric conformal prediction region
# estimated area
area.nonparametric <- function(region){
  if(class(region) != "list"){ 
    stop("Only appropriate for nonparametric conformal prediction region")
  }
  bins <- length(region); wn <- 1 / bins
  area <- 0
  for(i in 1:bins){
    nonpar.region <- region[[i]]
    area <- area + wn * as.numeric(rep(c(-1,1), 
      length(nonpar.region)/2) %*% nonpar.region)
  }
  area
}
area.nonparametric(nonparaCI)
# local coverage
local.coverage(region = nonparaCI, data = data.readme, d = p, 
  nonparametric = "TRUE", bins = bins, at.data = "TRUE")
# marginal coverage
local.coverage(region = nonparaCI, data = data.readme, d = p, 
  nonparametric = "TRUE", bins = 1, at.data = "TRUE")

## LSLW conformal prediction region
# estimated area
mean(apply(LSLW, 1, diff))
# local coverage
local.coverage(region = LSLW, data = data.readme, d = p, 
  bins = bins, at.data = "TRUE")
# marginal coverage
local.coverage(region = LSLW, data = data.readme, d = p, 
  bins = 1, at.data = "TRUE")

## HD region
# estimated area
mean(apply(HDCI, 1, diff))
# local coverage
local.coverage(region = HDCI, data = data.readme, d = p, 
  bins = bins, at.data = "TRUE")
# marginal coverage
local.coverage(region = HDCI, data = data.readme, d = p, 
  bins = 1, at.data = "TRUE")




#########################################
## make plot (.png)
png(filename="gammasimexample.png")
par(mfrow = c(2,2), oma = c(4,4,0,0), mar = c(1,1,1,1))

# parametric conformal prediction region
ix <- sort(x, index.return = TRUE)$ix
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y), max(y)))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], type = "l", col = "red")
lines(x[ix], paraCI[ix, 2], type = "l", col = "red")
axis(2)

# nonparametric conformal prediction region
plot.nonparametric <- function(region, x, y, bins){
  if(class(region) != "list"){ 
    stop("Only appropriate for nonparametric conformal prediction region")
  }
  plot.new()
  plot.window(xlim = c(0,1), ylim = c(min(y), max(y)))
  points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
  for(i in 1:bins){ 
    nonpar.i  <- nonparaCI[[i]]
    odd <- which(1:length(nonpar.i) %% 2 == 1) 
    even <- which(1:length(nonpar.i) %% 2 == 0)
    segments(x0 = 1/bins * (i-1), y0 = nonpar.i, x1 =1/bins * i, 
      col = "red")
    if(i != 1){ 
      nonpar.i.previous <- nonparaCI[[i-1]]
      nonpar <- sort(c(nonpar.i, nonpar.i.previous))
      odd2 <- which(1:length(nonpar) %% 2 == 1) 
      even2 <- which(1:length(nonpar) %% 2 == 0)
      segments(x0 = 1/bins * (i-1), y0 = nonpar[odd2], 
        y1 = nonpar[even2], col = "red")
    }
  }
}    
plot.nonparametric(nonparaCI, x = x, y = y, bins = bins)

# least squares conformal prediction region
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,max(y)))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLW[ix, 1], type = "l", col = "red")
lines(x[ix], LSLW[ix, 2], type = "l", col = "red")
axis(1); axis(2)

# highest density region
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,max(y)))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], HDCI[ix, 1], type = "l", col = "red")
lines(x[ix], HDCI[ix, 2], type = "l", col = "red")
axis(1)

# axis labels
mtext("x", side = 1, line = 2.5, outer = TRUE, cex = 1.25)
mtext("y", side = 2, line = 2.5, outer = TRUE, cex = 1.25)


dev.off()




#########################################
## make plot (.pdf)
pdf(file="gammasimexample.pdf")
par(mfrow = c(2,2), oma = c(4,4,0,0), mar = c(1,1,1,1))

# parametric conformal prediction region
ix <- sort(x, index.return = TRUE)$ix
plot.new()
plot.window(xlim = c(0,1), ylim = c(min(y), max(y)))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], paraCI[ix, 1], type = "l", col = "red")
lines(x[ix], paraCI[ix, 2], type = "l", col = "red")
axis(2)

# nonparametric conformal prediction region
plot.nonparametric <- function(region, x, y, bins){
  if(class(region) != "list"){ 
    stop("Only appropriate for nonparametric conformal prediction region")
  }
  plot.new()
  plot.window(xlim = c(0,1), ylim = c(min(y), max(y)))
  points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
  for(i in 1:bins){ 
    nonpar.i  <- nonparaCI[[i]]
    odd <- which(1:length(nonpar.i) %% 2 == 1) 
    even <- which(1:length(nonpar.i) %% 2 == 0)
    segments(x0 = 1/bins * (i-1), y0 = nonpar.i, x1 =1/bins * i, 
      col = "red")
    if(i != 1){ 
      nonpar.i.previous <- nonparaCI[[i-1]]
      nonpar <- sort(c(nonpar.i, nonpar.i.previous))
      odd2 <- which(1:length(nonpar) %% 2 == 1) 
      even2 <- which(1:length(nonpar) %% 2 == 0)
      segments(x0 = 1/bins * (i-1), y0 = nonpar[odd2], 
        y1 = nonpar[even2], col = "red")
    }
  }
}    
plot.nonparametric(nonparaCI, x = x, y = y, bins = bins)

# least squares conformal prediction region
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,max(y)))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], LSLW[ix, 1], type = "l", col = "red")
lines(x[ix], LSLW[ix, 2], type = "l", col = "red")
axis(1); axis(2)

# highest density region
plot.new()
plot.window(xlim = c(0,1), ylim = c(0,max(y)))
points(x, y, pch = 19, col = rgb(0,0,0,alpha=0.2))
lines(x[ix], HDCI[ix, 1], type = "l", col = "red")
lines(x[ix], HDCI[ix, 2], type = "l", col = "red")
axis(1)

# axis labels
mtext("x", side = 1, line = 2.5, outer = TRUE, cex = 1.25)
mtext("y", side = 2, line = 2.5, outer = TRUE, cex = 1.25)


dev.off()



file <- paste("readme", collapse = "")
file <- paste(file, "RData", sep = ".", collapse = "")
save.image(file = file, ascii = TRUE)

