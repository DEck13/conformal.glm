
find.index <- function (mat, wn, d) 
{
    n.out <- nrow(mat)
    indices <- rep(1, n.out)
    if (wn < 1) {
        A <- seq(from = 0, to = 1 - wn, by = wn)
        mat <- apply(mat, 2, function(x) {
            if (min(x) < 0 || max(x) > 1) {
                x <- (x - min(x)) + 1e-04
                x <- x/(sign(max(x)) * max(x)) 
            }
            x
        })
        for (j in 1:n.out) {
            foo <- 0
            for (i in 1:d) {
                if (i < d) {
                  foo <- foo + (max(which(A < mat[j, i])) - 1) * 
                    (1/wn)^(d - i)
                }
                if (i == d) {
                  foo <- foo + max(which(A < mat[j, i]))
                }
            }
            indices[j] <- foo
        }
    }
    indices
}








local.coverage <- function(region, nonparametric = "FALSE", data, 
  at.data = "TRUE", newdata = NULL, d, bins = NULL){

  n <- nrow(data)
  wn <- min(1/ floor(1 / (log(n)/n)^(1/(d+3))), 1/2)
  if(class(bins) != "NULL") wn <- 1 / bins
  index.data <- find.index(as.matrix(data[, 1:d + 1]), wn = wn, d = d)

  # parametric and nonparametric prediction regions are 
  # data is formatted as (y, predictors)
  output <- rep(NA, bins)
  if(at.data == TRUE){
    if(nonparametric == "FALSE"){
      x <- as.matrix(data[,1:d + 1], col = d)
      index.newdata <- find.index(x, wn = wn, d = d)
      colnames(region) <- c("lwr","upr")
      output[sort(unique(index.data))] <- 
        unlist(lapply(sort(unique(index.data)), function(j){    
          Y <- data[index.data == j, 1]
          lwr <- region[index.newdata == j, 1]
          upr <- region[index.newdata == j, 2]
          out <- mean(lwr <= Y & Y <= upr)
          out
      }))
    }
    if(nonparametric == "TRUE"){

      f <- sapply(data, is.factor)
      any.factors <- !all(f == FALSE)

      ## when there are no factors
      if(!any.factors){
        n.bins.region <- length(region)
        x <- as.matrix(data[,1:d + 1], col = d)
        index.bins.region <- find.index(x, wn = 1/n.bins.region, d = d)
        Y <- data[, 1]
        ## put each y value with its corresponding prediction interval
        foo <- lapply(1:n, FUN = function(j){
          c(Y[j], region[[index.bins.region[j]]])
        })

        output[sort(unique(index.data))] <- 
          unlist(lapply(sort(unique(index.data)), function(j){
            index.j <- which(index.data == j)
            ## indicate which endpoints y is between
            int <- unlist(lapply(index.j, FUN = function(x){
              bar <- 0
              if(any(foo[[x]][1] > foo[[x]][-1])){
                bar <- max(which(foo[[x]][1] > foo[[x]][-1]))
              }
              bar
            }))
            ## only odd numbers correspond to y being in one of the 
            ## (possibly disjoint) intervals
            mean(int %% 2 == 1)
        }))
      }

      ## when there are factors
      if(any.factors){
        n.factor.combinations <- length(region)
        d <- ncol(data) - 1
        ## data assumed to be of the form (response, predictors)
        X <- as.matrix(data[,1:d + 1], col = d)
        Y <- data[, 1]

        ## location of factor and numeric variables in dataframe 
        ## with response removed
        index.factor.variables <- which(f) - 1
        index.numeric.variables <- which(!f)[-1] - 1
        factors <- lapply(index.factor.variables, 
          function(j) as.numeric(as.factor(X[, j])))

        ## split model matrix by factor level combinations
        split.X.factors <- split(X, factors, drop = TRUE)

        ## split bin indices by factor level combinations
        bin.index.by.factors <- lapply(split.X.factors, function(x){ 
          mat <- matrix(x, ncol = ncol(X))
          colnames(mat) <- colnames(X)
          mat <- mat[, index.numeric.variables]
          find.index(mat = mat, wn = wn, d = ncol(mat))
        })

        ## split response by factor level combinations
        split.Y.factors <- split(Y, factors, drop = TRUE)
        resps.by.factors <- lapply(split.Y.factors, as.numeric)

        ## compute coverage probability
        output <- mean(unlist(lapply(1:n.factor.combinations, function(i){
          sapply(unique(bin.index.by.factors[[i]]), function(j){
            Y.int <- resps.by.factors[[i]][bin.index.by.factors[[i]] == j]
            region.int <- region[[i]][[j]]
            sapply(Y.int, function(y){
              inclusion <- 0
              if(any(y - region.int == 0)){
                inclusion <- 1
                return(inclusion)
              }
              inclusion.index <- y - region.int > 0
              if(!any(inclusion.index)){ 
                return(inclusion)
              }
              if(max(which(inclusion.index)) %% 2 == 1){ 
                inclusion <- 1
                return(inclusion)
              }
              return(inclusion)         
            })
          })
        })))

      }
    }
  }

  if(at.data == FALSE){
    ### right now only works for univariate regression
    ### and unimodal regions
    colnames(region) <- c("lwr","upr")
    index.newdata <- find.index(newdata, wn = wn, d = d)
    output[sort(unique(index.data))] <- 
      unlist(lapply(sort(unique(index.data)), function(j){
        input1 <- data[, 1:d + 1]
        m.lwr <- lm(lwr ~ poly(input1, degree = 3), 
          data = as.data.frame(cbind(data, region)), 
          subset = index.newdata == j)
        m.upr <- lm(upr ~ poly(input1, degree = 3), 
          data = as.data.frame(cbind(data, region)),
          subset = index.newdata == j)

        newdat <- as.data.frame(data[index.data == j, 1:d + 1])
        colnames(newdat) <- paste("input", 1:d, sep = "")
        p.lwr <- predict(m.lwr, newdata = newdat)
        p.upr <- predict(m.upr, newdata = newdat)
        y <- data[index.data == j, 1]
        out <- mean(p.lwr <= y & y <= p.upr)
    }))  
  }
  output
}




## conformal scores
phatxy <- function(ynew, xnew, Yk, Xk, xnew.modmat, 
  data, formula, family, link){

  out <- rateMLE.y <- shapeMLE.y <- NULL
  data.y <- rbind(data, data[1, ])
  n <- nrow(data)
  data.y[n+1,] <- cbind(ynew, xnew)

  if(family == "Gamma"){
    if(link == "identity"){
      m1.y <- glm(formula, data = data.y, family = Gamma(link = identity), 
        control = list(maxit = 1e4))
      shapeMLE.y <- as.numeric(gamma.shape(m1.y)[1])              
      rateMLE.y <- 1 / (cbind(1, rbind(Xk, xnew.modmat)) %*% 
        coefficients(m1.y)) * shapeMLE.y
    }
    if(link == "inverse"){
      m1.y <- glm(formula, data = data.y, family = Gamma(link = inverse), 
        control = list(maxit = 1e4))
      shapeMLE.y <- as.numeric(gamma.shape(m1.y)[1])              
      rateMLE.y <- (cbind(1, rbind(Xk, xnew.modmat)) %*% 
        coefficients(m1.y)) * shapeMLE.y
    }
    if(link == "log"){
      m1.y <- glm(formula, data = data.y, family = Gamma(link = log), 
        control = list(maxit = 1e4))
      shapeMLE.y <- as.numeric(gamma.shape(m1.y)[1])              
      rateMLE.y <- (1 / exp(cbind(1, rbind(Xk, xnew.modmat)) %*% 
        coefficients(m1.y))) * shapeMLE.y
    }

    out <- dgamma(c(Yk, ynew), rate = rateMLE.y, shape = shapeMLE.y)
  }

  if(family == "gaussian"){
    m1.y <- lm(formula, data = data.y)
    out <- dnorm(c(Yk, ynew), mean = as.numeric(cbind(1, rbind(Xk, xnew.modmat)) %*% 
      coefficients(m1.y)), sd = summary(m1.y)$sigma)
  }

  if(family == "inverse.gaussian"){
    m1.y <- glm(formula, data = data.y, family = family)
    out <- dinvgauss(c(Yk, ynew), mean = 1 / sqrt(cbind(1, rbind(Xk, xnew.modmat)) %*% 
      coefficients(m1.y)))
  }

  out
}






regions <- function(formula, data, newdata, family = "gaussian", link, 
  alpha = 0.10, cores = 1, bins = 1, parametric = TRUE, 
  nonparametric = FALSE, h = NULL, precision = 0.005){

  ## initial quantities
  respname <- all.vars(formula)[1]
  Y <- data[, colnames(data) %in% respname]
  n <- length(Y)
  n.pred <- nrow(newdata)

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
  X.variables <- as.matrix(data[,-which(colnames(data) %in% respname)], 
    nrow = n)
  d <- ncol(X.variables)
  m1 <- X <- betaOLS <- NULL; p <- sd.res <- 0
  if(family == "gaussian"){
  m1 <- lm(formula, data = data, x = TRUE)
    X <- matrix(m1$x[, -1], nrow = n)
    betaOLS <- betaMLE <- coefficients(m1)
    p <- length(betaOLS) - 1
    sd.res <- summary(m1)$sigma
  }
  ## Get MLEs and plugin interval for gamma distribution 
  betaMLE <- shapeMLE <- rateMLE <- 0
  if(family == "Gamma"){
    m1 <- glm(formula, data = data, family = family, x = TRUE)
    X <- matrix(m1$x[, -1], nrow = n)
    betaMLE <- coefficients(m1)
    p <- length(betaMLE) - 1
    shapeMLE <- as.numeric(gamma.shape(m1)[1])
    rateMLE <- cbind(1, X) %*% betaMLE * shapeMLE
  }

  ## Get MLEs and plugin interval for Inverse Gaussian Distribution 
  if(family == "inverse.gaussian"){
    m1 <- glm(formula, data = data, family = family, x = TRUE)
    X <- matrix(m1$x[, -1], nrow = n)
    betaMLE <- coefficients(m1)
    p <- length(betaMLE) - 1
  }

  ## set up partition
  wn <- min(1/ floor(1 / (log(n)/n)^(1/(d+3))), 1/2)
  if(class(bins) != "NULL") wn <- 1 / bins

  ## important internal quantities for our prediction regions   
  ## create seperate data frames with respect to each partition 
  ## ignores the intercept of the newdata matrix
  index <- find.index(X.variables, wn = wn, d = d)
  index.pred <- find.index(matrix(newdata, ncol = d), wn = wn, d = d)
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
        xnew <- matrix(newdata.variables[j, ], nrow = 1, ncol = d)
        xnew.modmat <- matrix(newdata.formula[j, ], nrow = 1, ncol = p)
        index.bin <- which(index == index.pred[j])
        nk <- length(index.bin)
        Xk <- matrix(X[index.bin, ], ncol = p)
        Yk <- Y[index.bin]

        ## conformal scores (density estimator)
        density_scores <- function(ynew){
          phatxy(ynew, xnew = xnew, Yk = Yk, Xk = Xk, xnew.modmat = xnew.modmat, 
            data = data, formula = formula, family = family, link = link)
        }

        ## declare output variable
        interval <- NULL
        
        ## initial check for enough data within bin
        nk.tilde <- floor(alpha * (nk + 1))
        if(nk.tilde <= 1){ 
           interval <- c(min(Yk), max(Yk))
        }

        ## set up a lower (upper lower bound) and upper bound 
        ## (lower upper bound) to start two line searchs in order to 
        ## construct the parametric conformal prediction region
        if(nk.tilde > 1){
          quant.Yk <- quantile(Yk, probs = c(alpha, 1 - alpha ))
          y.lwr <- y.min <- as.numeric(quant.Yk[1])
          y.upr <- y.max <- as.numeric(quant.Yk[2])

          ## lower line search
          prec <- max( min(diff(sort(Yk[Yk <= y.min]))), precision)      
          steps <- 1
          flag <- FALSE
          while(rank(density_scores(y.lwr))[nk + 1] >= nk.tilde & flag == FALSE){
            y.lwr <- y.lwr - steps * prec
            if(family != "gaussian"){ 
              if(y.lwr <= 0.0001){ 
                y.lwr <- 0.0001
                flag <- TRUE
              }
            }
            steps <- steps + 1
          }
          if(flag == FALSE){
            steps <- 1
            if(y.min - y.lwr <= 0.001){ 
              prec <- min( min(diff(sort(Yk[Yk <= y.min]))), precision) / 2
              if(prec == 0) prec <- precision / 2
            }
            if(y.min - y.lwr > 0.001){
              prec <- min( min(diff(sort(Yk[Yk <= y.min]))), precision)
              if(prec == 0) prec <- precision
            }
            while(rank(density_scores(y.lwr))[nk + 1] < nk.tilde & y.lwr < max(Yk)){
              y.lwr <- y.lwr + steps * prec
              steps <- steps + 1
            }
          }

          ## upper line search
          if(y.lwr >= y.upr) y.upr <- 2^sign(max(Yk)) * max(Yk)
          prec <- max( min(diff(sort(Yk[Yk >= y.max]))), precision)
          steps <- 1
          while(rank(density_scores(y.upr))[nk + 1] >= nk.tilde){
            y.upr <- y.upr + steps * prec
            steps <- steps + 1
          }
          steps <- 1
          if(y.upr - y.max <= precision){ 
            prec <- min( min(diff(sort(Yk[Yk >= y.max]))), precision) / 2
            if(prec == 0) prec <- precision / 2
          }
          if(y.upr - y.max > precision){
            prec <- min( min(diff(sort(Yk[Yk >= y.max]))), precision)
            if(prec == 0) prec <- precision
          }        
          while(rank(density_scores(y.upr))[nk + 1] < nk.tilde){
            y.upr <- y.upr - steps * prec
            if(family != "gaussian"){ 
              if(y.upr <= 0.0001){ 
                y.upr <- 0.0001
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
          interval <- c(y.lwr, y.upr)
        }
        interval  
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
    ## Nonparametric conformal prediction region 
    if(class(h) == "NULL") h <- wn
    COPS <- function(index){ 
      out <- mclapply(sort(unique(index)), mc.cores = cores, 
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
            round(sum(dnorm(Yknotj, mean = y, sd = h) 
              - dnorm(Yknotj, mean = Ykj, sd = h)), 6)
          })) > 0)
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
              c(endpts, y.seq[breaks[k]], y.seq[breaks[k]+1])
            if(k == length(breaks)) endpts <- 
              c(endpts, y.seq[breaks[k]], y.seq[breaks[k]+1], max(y.seq))
          }
        } 
        endpts
      })
      out
    }
    nonparaconformal <- COPS(index.pred)
  }

  out = list(paraconformal = paraconformal, 
    nonparaconformal = nonparaconformal)
  return(out)
}


