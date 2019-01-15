
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





local.coverage <- function(region, nonparametric = "FALSE", data, 
  at.data = "TRUE", newdata = NULL, k, bins = NULL){

  n <- nrow(data)
  wn <- min(1/ floor(1 / (log(n)/n)^(1/(k+3))), 1/2)
  if(class(bins) != "NULL") wn <- 1 / bins
  index.data <- find.index(as.matrix(data[, 1:k + 1]), wn = wn, k = k)

  # parametric and nonparametric prediction regions are 
  # data is formatted as (y, predictors)
  output <- NULL
  if(at.data == TRUE){
    if(nonparametric == "FALSE"){
      x <- as.matrix(data[,1:k + 1], col = k)
      index.newdata <- find.index(x, wn = wn, k = k)
      colnames(region) <- c("lwr","upr")
      output <- unlist(lapply(sort(unique(index.data)), function(j){    
        y <- data[index.data == j, 1]
        lwr <- region[index.newdata == j, 1]
        upr <- region[index.newdata == j, 2]
        out <- mean(lwr <= y & y <= upr)
        out
      }))
    }
    if(nonparametric == "TRUE"){
      n.bins.region <- length(region)
      x <- as.matrix(data[,1:k + 1], col = k)
      index.bins.region <- find.index(x, wn = 1/n.bins.region, k = k)
      y <- data[, 1]
      ## put each y value with its corresponding prediction interval
      foo <- lapply(1:n, FUN = function(j){
        c(y[j], region[[index.bins.region[j]]])
      })

      output <- unlist(lapply(sort(unique(index.data)), function(j){
        index.j <- which(index.data == j)
        ## indicate which endpoints y is between
        int <- unlist(lapply(index.j, FUN = function(x){
          which(foo[[x]][1] > foo[[x]][-1])
        }))
        ## only odd numbers correspond to y being in one of the 
        ## (possibly disjoint) intervals
        mean(int %% 2 == 1)
      }))
    }
  }

  if(at.data == FALSE){
    ### right now only works for univariate regression
    ### and unimodal regions
    colnames(region) <- c("lwr","upr")
    index.newdata <- find.index(newdata, wn = wn, k = k)
    output <- unlist(lapply(sort(unique(index.data)), function(j){
      input1 <- data[, 1:k + 1]
      m.lwr <- lm(lwr ~ poly(input1, degree = 3), 
        data = as.data.frame(cbind(data, region)), 
        subset = index.newdata == j)
      m.upr <- lm(upr ~ poly(input1, degree = 3), 
        data = as.data.frame(cbind(data, region)),
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
  nonparametric = FALSE, h = NULL, precision = 0.001){

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
        prec <- max( min(diff(sort(Yk[Yk <= y.min]))), precision)      
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
            prec <- min( min(diff(sort(Yk[Yk <= y.min]))), precision) / 2
          }
          if(y.min - y.lwr > 0.001){
            prec <- min( min(diff(sort(Yk[Yk <= y.min]))), precision)
            #if(prec < precision){ 
            #  prec <- mean(min(diff(sort(Yk[Yk <= y.min]))), 0.0001)
            #}
          }
          while(rank(phatxy(y.lwr))[nk + 1] < nk.tilde & y.lwr < max(Yk)){
            y.lwr <- y.lwr + steps * prec
            steps <- steps + 1
          }
        }

        ## upper line search
        if(y.lwr >= y.upr) y.upr <- 2^sign(max(Yk)) * max(Yk)
        prec <- max( min(diff(sort(Yk[Yk >= y.max]))), precision)
        steps <- 1
        while(rank(phatxy(y.upr))[nk + 1] >= nk.tilde){
          y.upr <- y.upr + steps * prec
          steps <- steps + 1
        }
        steps <- 1
        if(y.upr - y.max <= precision){ 
          prec <- min( min(diff(sort(Yk[Yk >= y.max]))), precision)
        }
        if(y.upr - y.max > precision){
          prec <- min( min(diff(sort(Yk[Yk >= y.max]))), precision)
          #if(prec < 0.001){ 
          #  prec <- mean(min(diff(sort(Yk[Yk >= y.max]))), 0.0001)
          #}
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
    if(class(h) == "NULL") h <- wn
    
    #COPS <- function(newdata){ 
    #  out <- mclapply(1:n.pred, mc.cores = cores, FUN = function(j){
    #    index.bin <- which(index == index.pred[j])
    #    nk <- length(index.bin)
    #    Yk <- Y[index.bin]

    #    ## initial check for enough data within bin
    #    nk.tilde <- floor(alpha * (nk + 1))
    #   if(nk.tilde == 0) stop("bin width is too small")

    #    ## nonparametric density
    #    phatxy <- function(y){                
    #      out <- which(unlist(lapply(1:nk, FUN = function(j){
    #       Yknotj <- Yk[-j]
    #        Ykj <- Yk[j]
    #        sum(dnorm(Yknotj, mean = y, sd = h) 
    #          - dnorm(Yknotj, mean = Ykj, sd = h))
    #      })) >= 0)
    #      if(length(out) == 0) out <- -1
    #      length(out) >= nk.tilde
    #    }

    #    ## set up a lower (upper lower bound) and upper bound 
    #    ## (lower upper bound) to start two line searchs in order to 
    #    ## construct the nonparametric conformal prediction region
    #    quant.Yk <- quantile(Yk, probs = c(2 * alpha, 1 - 2 * alpha ))
    #    y.lwr <- y.min <- as.numeric(quant.Yk[1])
    #    y.upr <- y.max <- as.numeric(quant.Yk[2])

    #    # lower line search
    #    prec <- max( min(diff(sort(Yk[Yk <= y.min]))), 0.001)      
    #    steps <- 1
    #    while(phatxy(y.lwr)){
    #      y.lwr <- y.lwr - steps * prec
    #      steps <- steps + 1
    #    }
    #    steps <- 1
    #    prec <- min( min(diff(sort(Yk[Yk <= y.min]))), 0.001)
    #    if(prec < 0.001) prec <- mean(min(diff(sort(Yk[Yk <= y.min]))), 0.001)
    #    while(!phatxy(y.lwr)){
    #      y.lwr <- y.lwr + steps * prec
    #      steps <- steps + 1
    #    }

    #    # upper line search
    #    prec <- max( min(diff(sort(Yk[Yk >= y.max]))), 0.001)
    #    steps <- 1
    #    while(phatxy(y.upr)){
    #      y.upr <- y.upr + steps * prec
    #      steps <- steps + 1
    #    }
    #    steps <- 1
    #    prec <- min( min(diff(sort(Yk[Yk >= y.max]))), 0.001) 
    #    if(prec < 0.001) prec <- mean(min(diff(sort(Yk[Yk >= y.max]))), 0.001)
    #    while(!phatxy(y.upr)){
    #      y.upr <- y.upr - steps * prec
    #      steps <- steps + 1
    #    }

    #    c(y.lwr, y.upr)
    #  })

    #  out <- do.call(rbind, out)
    #  colnames(out) <- c("lwr", "upr")
    #  out
    #}

    #nonparaconformal <- COPS(newdata)

    COPS <- function(newdata){ 
      out <- mclapply(sort(unique(index.pred)), mc.cores = cores, 
        FUN = function(j){

        index.bin <- which(index == j)
        nk <- length(index.bin)
        Yk <- Y[index.bin]

        ## initial check for enough data within bin
        nk.tilde <- floor(alpha * (nk + 1))
        if(nk.tilde == 0) stop("bin width is too small")

        ## nonparametric density
        phatxy <- function(y){                
          int <- which(unlist(lapply(1:nk, FUN = function(j){
            Yknotj <- Yk[-j]
            Ykj <- Yk[j]
            sum(dnorm(Yknotj, mean = y, sd = h) 
              - dnorm(Yknotj, mean = Ykj, sd = h))
          })) >= 0)
          if(length(int) == 0) int <- -1
          out <- length(int) >= nk.tilde
          out
        }

        # lower line search
        y.lwr <- y.min <- min(Yk)
        quant.Yk <- quantile(Yk, probs = c(2 * alpha, 1 - 2 * alpha ))
        prec <- max( min(diff(sort(Yk[Yk <= quant.Yk[1]]))), precision)      
        steps <- 1
        while(!phatxy(y.lwr)){
          y.lwr <- y.lwr + steps * prec
          steps <- steps + 1
        }
        steps <- 1
        prec <- min( min(diff(sort(Yk[Yk <= quant.Yk[1]]))), precision)
        #if(prec < 0.001) prec <- 0.0005
          #mean(min(diff(sort(Yk[Yk <= quant.Yk[1]]))), 0.001)
        while(phatxy(y.lwr)){
          y.lwr <- y.lwr - steps * prec
          steps <- steps + 1
        }

        # upper line search
        y.upr <- y.max <- max(Yk)
        prec <- max( min(diff(sort(Yk[Yk >= quant.Yk[2]]))), precision)
        steps <- 1
        while(!phatxy(y.upr)){
          y.upr <- y.upr - steps * prec
          steps <- steps + 1
        }
        steps <- 1
        prec <- min( min(diff(sort(Yk[Yk >= quant.Yk[2]]))), precision)
        #if(prec < 0.001) prec <- 0.0005
          #mean(min(diff(sort(Yk[Yk >= quant.Yk[2]]))), 0.001)
        while(phatxy(y.upr)){
          y.upr <- y.upr + steps * prec
          steps <- steps + 1
        }

        y.seq <- seq(y.lwr, y.upr, by = precision)
        foo <- unlist(lapply(y.seq, FUN = phatxy))
        y.seq <- y.seq[foo]
        breaks <- which(round(diff(y.seq), 
          ceiling(log10(1/precision))) != precision)
        endpts <- c(min(y.seq), max(y.seq))
        if(length(breaks) == 1) endpts <- 
          c(min(y.seq), y.seq[breaks], y.seq[breaks+1],  max(y.seq))
        if(length(breaks) > 1){
          for(k in 1:length(breaks)){
            if(k == 1) endpts <- 
              c(min(y.seq), y.seq[breaks[k]], y.seq[breaks[k]+1])
            if(k != 1 & k != length(breaks)) endpts <- 
              c(endpts, y.seq[breaks[k]+1], y.seq[breaks[k+1]])
            if(k == length(breaks)) endpts <- 
              c(endpts, y.seq[breaks[k]], y.seq[breaks[k]+1], max(y.seq))
          }
        } 
        endpts
      })
      out
    }
    nonparaconformal <- COPS(newdata)
  }

  out = list(paraconformal = paraconformal, 
    nonparaconformal = nonparaconformal)
  return(out)
}


