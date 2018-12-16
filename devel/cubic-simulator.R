

library(parallel)
library(MASS)
library(statmod)
library(conformal.glm)
library(conformalInference)
library(HDInterval)
set.seed(13)


absolute.error <- function(y, region, squared = FALSE){
    lwr <- region[, 1]
    upr <- region[, 2]
    index <- which(!(lwr <= y & y <= upr))
    out <- sum(unlist(lapply(index, function(j){
      foo <- NULL
      if(squared == FALSE) foo <- min(abs(y[j] - lwr[j]), abs(y[j] - upr[j]))
      if(squared == TRUE) foo <- (min(abs(y[j] - lwr[j]), abs(y[j] - upr[j])))^2
      foo
    }))) / length(y)
    out
}


cubic.simulator <- function(n, alpha = 0.10, beta, bins = NULL, family = "Gamma", 
  link = "inverse", shape = NULL, sd = NULL, confamily = "gaussian",
  parametric = TRUE, nonparametric = FALSE, 
  LS = FALSE, HD = FALSE, cores = 6){

	p <- k <- length(beta) - 1
  x <- matrix(runif(n), ncol = p)
  y <- rep(0,n)
  data <- NULL

  ## set up partition
  wn <- min(1/ floor(1 / (log(n)/n)^(1/(k+3))), 1/2)
  if(class(bins) == "NULL") bins <- 1 / wn


  if(family == "Gamma"){
    if(link == "identity"){
      rate <- (1 / (cbind(1, x) %*% beta)) * shape
      y <- rgamma(n = n, shape = shape, rate = rate)
      y <- y / sd(y)
      data <- data.frame(y = y, x = x)
      colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")
    }
    if(link == "inverse"){
      rate <- (cbind(1, x) %*% beta) * shape
      y <- rgamma(n = n, shape = shape, rate = rate)
      y <- y / sd(y)
      data <- data.frame(y = y, x = x)
      colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")
    }
    if(link == "log"){
      rate <- (1 / exp(cbind(1, x) %*% beta)) * shape
      y <- rgamma(n = n, shape = shape, rate = rate)
      y <- y / sd(y)
      data <- data.frame(y = y, x = x)
      colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")
    }
  }

  if(family == "gaussian"){
    mu <- cbind(1, x) %*% beta
    y <- rnorm(n = n, mean = mu, sd = sd)
    data <- data.frame(y = y, x = x)
    colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")
  }

  if(family == "inverse.gaussian"){
    mu = 1 / sqrt(cbind(1, x) %*% beta)
    y <- rinvgauss(n = n, mean = mu)
    y <- y / sd(y)
    data <- data.frame(y = y, x = x)
    colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")
  }

  tr <- try(fit <- glm(y ~ x1 + I(x1^2) + I(x1^3), family = "gaussian", data = data))
  paraCI <- nonparaCI <- LSCI <- HDCI <- NULL
  newdata <- NULL

  if(class(tr) != "try-error"){
    formula <- fit$formula
    newdata <- data
    respname <- all.vars(formula)[1]
    newdata <- newdata[, !(colnames(data) %in% respname)]
    newdata <- as.matrix(newdata)
  }
  if(class(tr) == "try-error"){
    parametric <- nonparametric <- LS <- HD <- FALSE
  }

  if(parametric){ 
    cpred <- conformal.glm(fit, parametric = TRUE, 
      nonparametric = FALSE, alpha = alpha,
      bins = bins, cores = cores)
    paraCI <- cpred$paraconformal
  }
  if(nonparametric){ 
    cpred <- conformal.glm(fit, parametric = FALSE, 
      nonparametric = TRUE, alpha = alpha,
      bins = bins, cores = cores)
    nonparaCI <- cpred$nonparaconformal
  }
  if(LS){
    funs <- lm.funs(intercept = TRUE)
    train.fun <- funs$train.fun
    predict.fun <- funs$predict.fun
    p1.tibs <- conformal.pred(x = cbind(x, x^2, x^3), y = y, x0 = cbind(x, x^2, x^3), 
      train.fun = train.fun, predict.fun = predict.fun, 
      alpha = alpha, grid.method = "linear",
      num.grid.pts = 999)
    LSCI <- cbind(p1.tibs$lo, p1.tibs$up)
  }
  if(HD){
    if(confamily == "Gamma"){
      betaMLE <- coefficients(fit)
      shapeMLE <- as.numeric(gamma.shape(fit)[1])
      rateMLE <- cbind(1, newdata) %*% betaMLE * shapeMLE
      HDCI <- do.call(rbind, lapply(1:nrow(newdata), function(j){ 
        hdi(qgamma, 1 - alpha, shape = shapeMLE, rate = rateMLE[j, 1])
      }))
    }
    if(confamily == "gaussian"){
      fit = lm(y ~ x1 + I(x1^2) + I(x1^3), data = data)
      betaMLE <- coefficients(fit)
      sdMLE <- summary(fit)$sigma
      meanMLE <- as.numeric(cbind(1, x, x^2, x^3) %*% betaMLE)
      HDCI <- do.call(rbind, lapply(1:nrow(newdata), function(j){ 
        hdi(qnorm, 1 - alpha, sd = sdMLE, mean = meanMLE[j])
      }))
    }
    if(confamily == "inverse.gaussian"){
      ## Blank for now
    }

  }

  ## local coverage and prediction error for prediction regions
  output.parametric <- output.nonparametric <- 
    output.LS <- output.HD <- rep(NA, bins + 1)
  if(parametric){
    local.parametric <- local.coverage(region = paraCI, 
      data = data, newdata = newdata, k = p, bins = bins, 
      at.data = "TRUE")
    para.pred.error <- absolute.error(y, paraCI, squared = TRUE)
    output.parametric <- c(local.parametric, para.pred.error, 
      mean(apply(paraCI, 1, diff)))
  }
  if(nonparametric){
    local.nonparametric <- local.coverage(region = nonparaCI, 
      data = data, newdata = newdata, k = p, bins = bins, 
      at.data = "TRUE")
    nonpara.pred.error <- absolute.error(y, nonparaCI, squared = TRUE)
    output.nonparametric <- c(local.nonparametric, nonpara.pred.error, 
      mean(apply(nonparaCI, 1, diff)))
  }
  if(LS){
    local.LS <- local.coverage(region = LSCI, 
      data = data, newdata = newdata, k = p, bins = bins, 
      at.data = "TRUE")
    LS.pred.error <- absolute.error(y, LSCI, squared = TRUE)
    output.LS <- c(local.LS, LS.pred.error, 
      mean(apply(LSCI, 1, diff))) 
  }
  if(HD){
    local.HD <- local.coverage(region = HDCI, 
      data = data, newdata = newdata, k = p, bins = bins, 
      at.data = "TRUE")
    HD.pred.error <- absolute.error(y, HDCI, squared = TRUE)
    output.HD <- c(local.HD, HD.pred.error, 
      mean(apply(HDCI, 1, diff))) 
  }

  output <- list(output.parametric = output.parametric, 
    output.nonparametric = output.nonparametric,
    output.LS = output.LS, 
    output.HD = output.HD)
  output

}





misgamgauss <- function(shape = shape, bins = bins){
  unlist(cubic.simulator(n = n, alpha = alpha, beta = beta, 
    bins = bins, family = "Gamma", 
    link = "inverse", shape = shape,
    confamily = "gaussian",
    parametric = TRUE, nonparametric = TRUE, 
    LS = TRUE, HD = TRUE, cores = 6))
}



## inputs
alpha <- 0.10
beta <- c(1, 5)
sd <- 2



## n = 500, bins = 2, shape = 1
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.1 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 1, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.1[, 1])))
B <- B - NAs
out500.2.1 <- cbind(apply(gamgauss500.2.1, 2, mean), 
  apply(gamgauss500.2.1, 2, sd) / sqrt(B))
rownames(out500.2.1) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.1


## n = 500, bins = 2, shape = 2
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.2 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 2, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.2[, 1])))
B <- B - NAs
out500.2.2 <- cbind(apply(gamgauss500.2.2, 2, mean), 
  apply(gamgauss500.2.2, 2, sd) / sqrt(B))
rownames(out500.2.2) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.2


## n = 500, bins = 2, shape = 5
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.5 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 5, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.5[, 1])))
B <- B - NAs
out500.2.5 <- cbind(apply(gamgauss500.2.5, 2, mean), 
  apply(gamgauss500.2.5, 2, sd) / sqrt(B))
rownames(out500.2.5) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.5


## n = 500, bins = 2, shape = 10
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.10 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 10, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.10[, 1])))
B <- B - NAs
out500.2.10 <- cbind(apply(gamgauss500.2.10, 2, mean), 
  apply(gamgauss500.2.10, 2, sd) / sqrt(B))
rownames(out500.2.10) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.10


## n = 500, bins = 2, shape = 25
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.25 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 25, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.25[, 1])))
B <- B - NAs
out500.2.25 <- cbind(apply(gamgauss500.2.25, 2, mean), 
  apply(gamgauss500.2.25, 2, sd) / sqrt(B))
rownames(out500.2.25) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.25


## n = 500, bins = 2, shape = 50
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.50 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 50, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.50[, 1])))
B <- B - NAs
out500.2.50 <- cbind(apply(gamgauss500.2.50, 2, mean), 
  apply(gamgauss500.2.50, 2, sd) / sqrt(B))
rownames(out500.2.50) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.50


## n = 500, bins = 2, shape = 100
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.100 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 100, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.100[, 1])))
B <- B - NAs
out500.2.100 <- cbind(apply(gamgauss500.2.100, 2, mean), 
  apply(gamgauss500.2.100, 2, sd) / sqrt(B))
rownames(out500.2.100) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.100


## n = 500, bins = 2, shape = 500
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.500 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 500, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.500[, 1])))
B <- B - NAs
out500.2.500 <- cbind(apply(gamgauss500.2.500, 2, mean), 
  apply(gamgauss500.2.500, 2, sd) / sqrt(B))
rownames(out500.2.500) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.500


## n = 500, bins = 2, shape = 1000
B <- 50
n <- 500
bins <- 2
system.time(gamgauss500.2.1000 <- do.call(rbind, lapply(1:B, 
  FUN = function(j){
    misgamgauss(shape = 1000, bins = bins)
  })))
NAs <- length(which(is.na(gamgauss500.2.1000[, 1])))
B <- B - NAs
out500.2.1000 <- cbind(apply(gamgauss500.2.1000, 2, mean), 
  apply(gamgauss500.2.1000, 2, sd) / sqrt(B))
rownames(out500.2.1000) <- c( paste("para.loc", 1:bins, sep = "."),
  "para.pred.error", "para.area",
  paste("nonpara.loc", 1:bins, sep = "."), 
  "nonpara.pred.error", "nonpara.area",
  paste("LS.loc", 1:bins, sep = "."), 
  "LS.pred.error", "LS.area", 
  paste("HD.loc", 1:bins, sep = "."), 
  "HD.pred.error", "HD.area")
out500.2.1000


file <- paste("cubic-simulator-output", collapse = "")
file <- paste(file, "RData", sep = ".", collapse = "")
save.image(file = file, ascii = TRUE)




## at data
para.area <- nonpara.area <- LS.area <- HD.area <- NULL
para.pred.error <- nonpara.pred.error <- LS.pred.error <- HD.pred.error <- NULL
para.local.coverage <- nonpara.local.coverage <- LS.local.coverage <- HD.local.coverage <- NULL
shapes <- c(1, 2, 5, 10, 25, 50, 100, 500, 1000)
shapes <- c(1, 2, 5, 10, 25, 50, 100, 500)
shapes <- c(1, 2, 5, 10, 25, 50, 100)
for(j in shapes ){
  foo <- eval(parse(text=paste("out500.2", j, sep = ".")))
  k <- which(shapes == j)
  para.local.coverage[c(2*k-1,2*k)] <- as.numeric(foo[1:2, 1])
  para.pred.error[k] <- as.numeric(foo[3, 1]) #* n
  para.area[k] <- as.numeric(foo[4, 1])
  nonpara.local.coverage[c(2*k-1,2*k)] <- as.numeric(foo[5:6, 1])
  nonpara.pred.error[k] <- as.numeric(foo[7, 1]) #* n
  nonpara.area[k] <- as.numeric(foo[8, 1])
  LS.local.coverage[c(2*k-1,2*k)] <- as.numeric(foo[9:10, 1])
  LS.pred.error[k] <- as.numeric(foo[11, 1]) #* n
  LS.area[k] <- as.numeric(foo[12, 1])
  HD.local.coverage[c(2*k-1,2*k)] <- as.numeric(foo[13:14, 1])
  HD.pred.error[k] <- as.numeric(foo[15, 1]) #* n
  HD.area[k] <- as.numeric(foo[16, 1])
}





png(filename="misspecif.png")
par(mfrow = c(2,2), oma = c(4,4,0,0), mar = c(1,1,1,1))


plot.new()
plot.window(xlim = c(log(min(shapes)), log(max(shapes))), ylim = c(0, 0.3))
lines(log(shapes), para.pred.error, col = "red", lty = 1)
lines(log(shapes), nonpara.pred.error, col = "green", lty = 2)
lines(log(shapes), LS.pred.error, col = "black", lty = 3)
lines(log(shapes), HD.pred.error, col = "blue", lty = 4)
axis(2)

legend(0.5, 0.30, legend=c("parametric conformal", "nonparametric conformal", "LS conformal", "HD region"),
       col=c("red", "green", "black", "blue"), lty=1:4)


plot.new()
plot.window(xlim = c(log(min(shapes)), log(max(shapes))), ylim = c(0, 3))
lines(log(shapes), para.area, col = "red", lty = 1)
lines(log(shapes), nonpara.area, col = "green", lty = 2)
lines(log(shapes), LS.area, col = "black", lty = 3)
lines(log(shapes), HD.area, col = "blue", lty = 4)
axis(2)

s <- length(shapes)
plot.new()
plot.window(xlim = c(log(min(shapes)), log(max(shapes))), ylim = c(0.80, 0.99))
lines(log(shapes), para.local.coverage[2*(1:s) - 1], col = "red", lty = 1)
lines(log(shapes), nonpara.local.coverage[2*(1:s) - 1], col = "green", lty = 2)
lines(log(shapes), LS.local.coverage[2*(1:s) - 1], col = "black", lty = 3)
lines(log(shapes), HD.local.coverage[2*(1:s) - 1], col = "blue", lty = 4)
axis(1); axis(2)

plot.new()
plot.window(xlim = c(log(min(shapes)), log(max(shapes))), ylim = c(0.85, 0.99))
lines(log(shapes), para.local.coverage[2*(1:s)], col = "red", lty = 1)
lines(log(shapes), nonpara.local.coverage[2*(1:s)], col = "green", lty = 2)
lines(log(shapes), LS.local.coverage[2*(1:s)], col = "black", lty = 3)
lines(log(shapes), HD.local.coverage[2*(1:s)], col = "blue", lty = 4)
axis(1); axis(2)

#dev.copy(png,'misspecif.png')
dev.off()

