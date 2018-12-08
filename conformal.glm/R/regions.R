
find.index <- function(mat, wn, k){
  A <- seq(from = 0, to = 1 - wn, by = wn)
  if(wn == 1) A <- 0
  n.out <- nrow(mat)
  out <- rep(0, n.out)
  for(j in 1:n.out){
    foo <- 0
    for(i in 1:k){
      if(i < k){
        foo <- foo + 
          (max(which(A < mat[j, i])) - 1) * (1/wn)^(k-i)
      }
      if(i == k){
        foo <- foo + 
          max(which(A < mat[j, i]))
      }
    }
    out[j] <- foo
  }
  out
}





local.coverage <- function(region, data, newdata, k, bins = NULL,
  at.data = "TRUE"){

  region <- cbind(newdata, region)
  n <- nrow(data)
  wn <- min(1/ floor(1 / (log(n)/n)^(1/(k+3))), 1/2)
  if(class(bins) != "NULL") wn <- 1 / bins

  index.data <- find.index(as.matrix(data[, 1:k + 1]), wn = wn, k = k)
  index.newdata <- find.index(newdata, wn = wn, k = k)
  # parametric and nonparametric prediction regions are 
  # formatted as (newdata, region).
  # data is formatted as (y, predictors)
  output <- NULL
  if(at.data == TRUE){
    colnames(region)[c(k+1,k+2)] <- c("lwr","upr")
    output <- unlist(lapply(sort(unique(index.data)), function(j){    
      y <- data[index.data == j, 1]
      lwr <- region[index.newdata == j, k+1]
      upr <- region[index.newdata == j, k+2]
      out <- mean(lwr <= y & y <= upr)
      out
    }))
  }

  if(at.data == FALSE){
    ### right now only works for univariate regression
    colnames(region)[c(k+1,k+2)] <- c("lwr","upr")
    output <- unlist(lapply(sort(unique(index.data)), function(j){
      input1 <- region[, 1:k]
      m.lwr <- lm(lwr ~ poly(input1, degree = 3), 
        data = as.data.frame(region)[, c(1:k, k+1)], 
        subset = index.newdata == j)
      m.upr <- lm(upr ~ poly(input1, degree = 3), 
        data = as.data.frame(region)[, c(1:k, k+2)], 
        subset = index.newdata == j)

      newdat <- as.data.frame(data[index.data == j, 1:k + 1])
      colnames(newdat) <- paste("input", 1:k, sep = "")
      p.lwr <- predict(m.lwr, newdata = newdat)
      p.upr <- predict(m.upr, newdata = newdat)
      y <- data[index.data == j, 1]
      out <- mean(p.lwr <= y & y <= p.upr)
    }))  
  }
  output
}





regions <- function(formula, data, newdata, family = "gaussian", link, 
  alpha = 0.10, cores = 1, bins = NULL, parametric = TRUE, 
  nonparametric = FALSE){

  ## initial quantities
  respname <- all.vars(formula)[1]
  Y <- data[, colnames(data) %in% respname]
  n <- length(Y)
  n.pred <- nrow(newdata)
  convergence <- 0

  if(is.null(newdata)){ 
    newdata <- data
    newdata <- newdata[, !(colnames(data) %in% respname)]
  }

  newdata <- as.matrix(newdata)
  colnames(newdata) <- colnames(data)[!(colnames(data) %in% respname)]
  ## create model calls when appropriate
  ## obtain OLS estimate 
  ## calculate important quantities for the gaussian 
  ## distribution 
  m1 <- lm(formula, data = data, x = TRUE)
  X <- matrix(m1$x[, -1], nrow = n)
  X.variables <- as.matrix(data[,-which(colnames(data) %in% respname)], 
    nrow = n)
  k <- ncol(X.variables)
  betaOLS <- betaMLE <- coefficients(m1)
  p <- length(betaOLS) - 1
  sd.res <- summary(m1)$sigma

  ## Get MLEs and plugin interval for gamma gamma distribution 
  shapeMLE <- rateMLE <- 0
  if(family == "Gamma"){
    m1 <- glm(formula, data = data, family = family)
    betaMLE <- coefficients(m1)
    shapeMLE <- as.numeric(gamma.shape(m1)[1])
    rateMLE <- cbind(1, X) %*% betaMLE * shapeMLE
  }

  ## Get MLEs and plugin interval for Inverse Gaussian Distribution 
  if(family == "inverse.gaussian"){
    m1 <- glm(formula, data = data, family = family)
    betaMLE <- coefficients(m1)
  }

  ## set up partition
  wn <- min(1/ floor(1 / (log(n)/n)^(1/(k+3))), 1/2)
  if(class(bins) != "NULL") wn <- 1 / bins

  ## important internal quantities for our prediction regions   
  ## create seperate data frames with respect to each partition 
  ## ignores the intercept of the newdata matrix
  index <- find.index(X.variables, wn = wn, k = k)
  index.pred <- find.index(matrix(newdata, ncol = k), wn = wn, k = k)
  indices.pred <- sort(unique(index.pred))


  ## newdata object for main effects only
  newdata.variables <- as.matrix(model.frame(~ ., as.data.frame(newdata)))

  ## newdata object for complete formula
  foo <- cbind(Y[1:nrow(newdata)], newdata)
  colnames(foo)[1] <- respname
  newdata.formula <- as.matrix(model.frame(formula, as.data.frame(foo))[, -1])
  

  paraconformal <- nonparaconformal <- NULL
  if(parametric == TRUE){

    # upfront quantities to improve speed
    shapeMLE.y <- shapeMLE; rateMLE.y <- rbind(rateMLE, 1); betaMLE.y <- betaMLE
    sd.y <- sd.res
    m1.y <- m1

    # ---- The parametric conformal implementation -------
    ## Parametric conformal prediction region for each 
    ## desired predictor combination 
    ## ignores the intercept of the newdata matrix
    Copt <- function(newdata){

      ## initial quantities
      data.y <- matrix(0, nrow = nrow(data) + 1, ncol = ncol(newdata.variables) + 1)
      colnames(data.y) <- colnames(data)

      out <- mclapply(1:n.pred, mc.cores = cores, FUN = function(j){
        x.variables <- matrix(newdata.variables[j, ], nrow = 1, ncol = k)
        x <- matrix(newdata.formula[j, ], nrow = 1, ncol = p)
        index.bin <- which(index == index.pred[j])
        nk <- length(index.bin)
        Xk <- matrix(X[index.bin, ], ncol = p)
        Yk <- Y[index.bin]

        ## conformal scores
        phatxy <- function(z){
          out <- rateMLE.y <- NULL
          data.y[, colnames(data) %in% respname] <- c(Y, z)
          data.y[, !(colnames(data) %in% respname)] <- rbind(X.variables, x.variables)
          data.y <- as.data.frame(data.y)

          if(family == "Gamma"){
            m1.y <- glm(formula, data = data.y, family = family)
            shapeMLE.y <- as.numeric(gamma.shape(m1.y)[1])

            if(link == "identity"){
              rateMLE.y <- 1 / (cbind(1, rbind(Xk, x)) %*% 
                coefficients(m1.y)) * shapeMLE.y
            }
            if(link == "inverse"){
              rateMLE.y <- (cbind(1, rbind(Xk, x)) %*% 
                coefficients(m1.y)) * shapeMLE.y
            }
            if(link == "log"){
              rateMLE.y <- (1 / exp(cbind(1, rbind(Xk, x)) %*% 
                coefficients(m1.y))) * shapeMLE.y
            }

            out <- dgamma(c(Yk, z), rate = rateMLE.y, shape = shapeMLE.y)
          }

          if(family == "gaussian"){
            m1.y <- lm(formula, data = data.y)
            out <- dnorm(c(Yk, z), mean = as.numeric(cbind(1, rbind(Xk, x)) %*% 
              coefficients(m1.y)), sd = summary(m1.y)$sigma)
          }

          if(family == "inverse.gaussian"){
            m1.y <- glm(formula, data = data.y, family = family)
            out <- dinvgauss(c(Yk, z), mean = 1 / sqrt(cbind(1, rbind(Xk, x)) %*% 
              coefficients(m1.y)))
          }

          out
        }

        ## initial check for enough data within bin
        nk.tilde <- floor(alpha * (nk + 1))
        if(nk.tilde == 0) stop("bin width is too small")

        ## set up a lower (upper lower bound) and upper bound 
        ## (lower upper bound) to start two line searchs in order to 
        ## construct the parametric conformal prediction region
        quant.Yk <- quantile(Yk, probs = c(alpha, 1 - alpha ))
        y.lwr <- y.min <- as.numeric(quant.Yk[1])
        y.upr <- y.max <- as.numeric(quant.Yk[2])


        ## lower line search
        prec <- max( min(diff(sort(Yk[Yk <= y.min]))), 0.001)      
        steps <- 1
        flag <- FALSE
        while(rank(phatxy(y.lwr))[nk + 1] >= nk.tilde & flag == FALSE){
          y.lwr <- y.lwr - steps * prec
          if(family != "gaussian"){ 
            if(y.lwr <= 0.00001){ 
              y.lwr <- 0.00001
              flag <- TRUE
            }
          }
          steps <- steps + 1
        }
        if(flag == FALSE){
          steps <- 1
          if(y.min - y.lwr <= 0.001){ 
            prec <- min( min(diff(sort(Yk[Yk <= y.min]))), 0.001) / 2
          }
          if(y.min - y.lwr > 0.001){
            prec <- min( min(diff(sort(Yk[Yk <= y.min]))), 0.001)
            if(prec < 0.001){ 
              prec <- mean(min(diff(sort(Yk[Yk <= y.min]))), 0.0001)
            }
          }
          while(rank(phatxy(y.lwr))[nk + 1] < nk.tilde & y.lwr < max(Yk)){
            y.lwr <- y.lwr + steps * prec
            steps <- steps + 1
          }
        }

        ## upper line search
        if(y.lwr >= y.upr) y.upr <- 2^sign(max(Yk)) * max(Yk)
        prec <- max( min(diff(sort(Yk[Yk >= y.max]))), 0.001)
        steps <- 1
        while(rank(phatxy(y.upr))[nk + 1] >= nk.tilde){
          y.upr <- y.upr + steps * prec
          steps <- steps + 1
        }
        steps <- 1
        if(y.upr - y.max <= 0.001){ 
          prec <- min( min(diff(sort(Yk[Yk >= y.max]))), 0.001)
        }
        if(y.upr - y.max > 0.001){
          prec <- min( min(diff(sort(Yk[Yk >= y.max]))), 0.001)
          if(prec < 0.001){ 
            prec <- mean(min(diff(sort(Yk[Yk >= y.max]))), 0.0001)
          }
        }        
        while(rank(phatxy(y.upr))[nk + 1] < nk.tilde){
          y.upr <- y.upr - steps * prec
          if(family != "gaussian"){ 
            if(y.upr < 0.00001){ 
              y.upr <- 0.00001
              break
            }
          }
          steps <- steps + 1
        }          
        if(abs(y.lwr - y.upr) < 0.002) y.lwr <- 0.00001
        if(y.lwr > y.upr){
          quant.Yk2 <- quantile(Yk, probs = c(alpha/2, 1 - alpha/2))
          y.lwr <- as.numeric(quant.Yk2[1])
          y.upr <- as.numeric(quant.Yk2[2])
        }

        c(y.lwr, y.upr)
      })

      #out <- t(out)
      out <- do.call(rbind, out)
      colnames(out) <- c("lwr", "upr")
      out
    }

    paraconformal <- Copt(newdata.variables)

  }


  if(nonparametric == TRUE){

    # ---- The nonparametric conformal implementation -------
    ## Nonparametric conformal prediction region for each 
    ## desired predictor combination 
    ## ignores the intercept of the newdata matrix
    #wn <- min(1/ floor(1 / (log(n)/n)^(1/(k+3))), 1/2)
    #hn <- ((log(n)/n)^(1/(1*(p+2)+1)))
    hn <- wn
    COPS <- function(newdata){ 
      out <- mclapply(1:n.pred, mc.cores = cores, FUN = function(j){
      #out <- apply(matrix(1:n.pred), 1, FUN = function(j){
        x.variables <- matrix(newdata.variables[j, ], nrow = 1, ncol = k)
        x <- matrix(newdata.formula[j, ], nrow = 1, ncol = p)
        index.bin <- which(index == index.pred[j])
        nk <- length(index.bin)
        Xk <- matrix(X[index.bin, ], ncol = p)
        Yk <- Y[index.bin]

        ## initial check for enough data within bin
        nk.tilde <- floor(alpha * (nk + 1))
        if(nk.tilde == 0) stop("bin width is too small")


        ## nonparametric density
        phatxy <- function(y){                
          out <- which(unlist(lapply(1:nk, FUN = function(j){
            Yknotj <- Yk[-j]
            Ykj <- Yk[j]
            sum(dnorm(Yknotj, mean = y, sd = hn) 
              - dnorm(Yknotj, mean = Ykj, sd = hn))
          })) >= 0)
          if(length(out) == 0) out <- -1
          length(out) >= nk.tilde
        }

        ## set up a lower (upper lower bound) and upper bound 
        ## (lower upper bound) to start two line searchs in order to 
        ## construct the parametric conformal prediction region
        quant.Yk <- quantile(Yk, probs = c(2 * alpha, 1 - 2 * alpha ))
        y.lwr <- y.min <- as.numeric(quant.Yk[1])
        y.upr <- y.max <- as.numeric(quant.Yk[2])

        # lower line search
        prec <- max( min(diff(sort(Yk[Yk <= y.min]))), 0.001)      
        steps <- 1
        while(phatxy(y.lwr)){
          y.lwr <- y.lwr - steps * prec
          steps <- steps + 1
        }
        steps <- 1
        prec <- min( min(diff(sort(Yk[Yk <= y.min]))), 0.001)
        if(prec < 0.001) prec <- mean(min(diff(sort(Yk[Yk <= y.min]))), 0.001)
        while(!phatxy(y.lwr)){
          y.lwr <- y.lwr + steps * prec
          steps <- steps + 1
        }

        # upper line search
        prec <- max( min(diff(sort(Yk[Yk >= y.max]))), 0.001)
        steps <- 1
        while(phatxy(y.upr)){
          y.upr <- y.upr + steps * prec
          steps <- steps + 1
        }
        steps <- 1
        prec <- min( min(diff(sort(Yk[Yk >= y.max]))), 0.001) 
        if(prec < 0.001) prec <- mean(min(diff(sort(Yk[Yk >= y.max]))), 0.001)
        while(!phatxy(y.upr)){
          y.upr <- y.upr - steps * prec
          steps <- steps + 1
        }

        c(y.lwr, y.upr)
      })

      out <- do.call(rbind, out)
      colnames(out) <- c("lwr", "upr")
      out
    }

    nonparaconformal <- COPS(newdata)
  }

  out = list(paraconformal = paraconformal, 
    nonparaconformal = nonparaconformal)
  return(out)
}


