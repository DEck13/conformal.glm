

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


simulator <- function(n, alpha = 0.10, beta, bins = NULL, family = "Gamma", 
  link = "inverse", shape = NULL, sd = NULL, confamily = "Gamma",
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
      data <- data.frame(y = y, x = x)
      colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")
    }
    if(link == "inverse"){
      rate <- (cbind(1, x) %*% beta) * shape
      y <- rgamma(n = n, shape = shape, rate = rate)
      data <- data.frame(y = y, x = x)
      colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")
    }
    if(link == "log"){
      rate <- (1 / exp(cbind(1, x) %*% beta)) * shape
      y <- rgamma(n = n, shape = shape, rate = rate)
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
    data <- data.frame(y = y, x = x)
    colnames(data)[2:(p+1)] <- paste("x", 1:p, sep = "")
  }

  tr <- try(fit <- glm(y ~ ., family = "gaussian", data = data))
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
    p1.tibs <- conformal.pred(x = x, y = y, x0 = newdata, 
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
      fit = lm(y ~ ., data = data)
      betaMLE <- coefficients(fit)
      sdMLE <- summary(fit)$sigma
      meanMLE <- as.numeric(cbind(1, newdata) %*% betaMLE)
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
    para.pred.error <- absolute.error(y, paraCI)
    output.parametric <- c(local.parametric, para.pred.error, 
      mean(apply(paraCI, 1, diff)))
  }
  if(nonparametric){
    local.nonparametric <- local.coverage(region = nonparaCI, 
      data = data, newdata = newdata, k = p, bins = bins, 
      at.data = "TRUE")
    nonpara.pred.error <- absolute.error(y, nonparaCI)
    output.nonparametric <- c(local.nonparametric, nonpara.pred.error, 
      mean(apply(nonparaCI, 1, diff)))
  }
  if(LS){
    local.LS <- local.coverage(region = LSCI, 
      data = data, newdata = newdata, k = p, bins = bins, 
      at.data = "TRUE")
    LS.pred.error <- absolute.error(y, LSCI)
    output.LS <- c(local.LS, LS.pred.error, 
      mean(apply(LSCI, 1, diff))) 
  }
  if(HD){
    local.HD <- local.coverage(region = HDCI, 
      data = data, newdata = newdata, k = p, bins = bins, 
      at.data = "TRUE")
    HD.pred.error <- absolute.error(y, HDCI)
    output.HD <- c(local.HD, HD.pred.error, 
      mean(apply(HDCI, 1, diff))) 
  }

  output <- list(output.parametric = output.parametric, 
    output.nonparametric = output.nonparametric,
    output.LS = output.LS, 
    output.HD = output.HD)
  output

}





