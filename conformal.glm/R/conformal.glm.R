

conformal.glm <- function(object, ..., newdata = NULL, alpha = 0.10, 
	cores = 1, bins = 1, parametric = TRUE, intercept = TRUE, 
	nonparametric = FALSE, h = NULL, precision = 0.005){

  out <- NULL 
  if(!any(attr(object$terms, "dataClasses")[-1] == "factor")){

    ## some important quantities
    call <- object$call
    formula <- call$formula
    fam <- object$family
    family <- fam$family
    link <- fam$link
    data <- object$data
    if(is.null(newdata)){ 
      newdata <- data
      respname <- all.vars(formula)[1]
      newdata <- newdata[, !(colnames(data) %in% respname)]
    }
    newdata <- as.matrix(newdata)

    stopifnot("glm" %in% class(object))
    stopifnot(family %in% c("Gamma", "gaussian", "inverse.gaussian"))

    int <- regions(formula = formula, data = data, 
      newdata = newdata, family = family, link = link, 
      alpha = alpha, cores = cores, bins = bins, 
      parametric = parametric, nonparametric = nonparametric, 
      h = h, precision = precision)

    paraconformal <- int$paraconformal
    nonparaconformal <- int$nonparaconformal
    out = list(paraconformal = paraconformal, 
      nonparaconformal = nonparaconformal)    
  }


  if(any(attr(object$terms, "dataClasses")[-1] == "factor")){

    object <- update(object, x = TRUE, y = TRUE)
    X <- object$x[, -1]
    Y <- object$y
    n <- length(Y)
    p <- ncol(X) - 1

    ## some important quantities
    call <- object$call
    formula <- call$formula
    respname <- all.vars(formula)[1]
    fam <- object$family
    family <- fam$family
    link <- fam$link
    data <- object$data
    terms <- c(respname, attr(terms(object), "term.labels"))
    data <- data[, colnames(data) %in% terms]
    d <- ncol(data)
    wn <- min(1/ floor(1 / (log(n)/n)^(1/(d+3))), 1/2)
    if(class(bins) != "NULL") wn <- 1 / bins

    stopifnot("glm" %in% class(object))
    stopifnot(family %in% c("Gamma", "gaussian", "inverse.gaussian"))

    newdata <- NULL
    if(is.null(newdata)){ 
      newdata <- data[, colnames(data) %in% attr(terms(object), "term.labels")]
    }
    n.pred <- nrow(newdata)

    ## model matrix for newdata
    maineffects.newdata <- model.matrix(~ ., newdata)[, -1]

    ## model matrix for newdata
    tt <- terms(object)
    Terms <- delete.response(tt)
    X.newdata <- model.matrix(~ ., model.frame(Terms, newdata))[, -1]

    ## split newdata model matrix by factor level combinations
    ## newdata is exactly the same as data but the first column is removed
    index.factor.variables <- which(attr(terms(object), "dataClasses") == "factor") - 1
    index.numeric.variables <- which(attr(terms(object), "dataClasses") == "numeric")[-1] - 1
    index.newdata.factor.variables <- which(colnames(newdata) %in% names(index.factor.variables))
    index.newdata.numeric.variables <- which(colnames(newdata) %in% names(index.numeric.variables))
    factors.newdata <- lapply(index.newdata.factor.variables, function(j) as.numeric(as.factor(newdata[, j])))
    split.newdata.factors <- split(newdata, factors.newdata, drop = FALSE)
    #preds.newdata.by.factors <- lapply(split.X.newdata.factors, function(x){ 
    #  mat <- matrix(x, ncol = ncol(X))
    #  colnames(mat) <- colnames(X)
    #  mat
    #})

    ## split newdata frame by factor level combinations
    newdata.by.factors <- split(newdata, factors.newdata, drop = FALSE)

    ## split newdata bin indices by factor level combinations
    newdata.bin.index.by.factors <- lapply(newdata.by.factors, function(x){ 
      if(nrow(x) == 0) return(0)
      mat <- as.matrix(x[, index.newdata.numeric.variables], ncol = length(index.newdata.numeric.variables))
      #colnames(mat) <- colnames(newdata[index.newdata.numeric.variables])
      output <- find.index(mat = mat, wn = wn, d = ncol(mat))
      output
    })

    ## split newdata model matrix by factor level combinations 
    ## and isolate numeric variables
    factor.code <- unlist(lapply(names(split.newdata.factors), function(x){
      paste(str_extract_all(x, "\\(?[0-9,]+\\)?")[[1]], collapse = "")    
    }))
    index.newdata <- unlist(lapply(1:nrow(newdata), function(x){ 
      index <- paste(as.numeric(newdata[x, index.newdata.factor.variables]), collapse = "")
      which(factor.code == index)
    }))

    ## need variable to say which bin the new data point belongs to
    index.bin.newdata <- unsplit(newdata.bin.index.by.factors, 
      factors.newdata, drop = FALSE)

    ## split response by factor level combinations
    #split.Y.factors <- split(Y, factors, drop = FALSE)
    #resps.by.factors <- lapply(split.Y.factors, as.numeric)

    ## get binning structure for newdata
    #nks.newdata <- as.numeric(lengths(numeric.newdata.by.factors)) / 
    #  length(index.numeric.variables)

    ## indices of nonempty factor level combinations in 
    ## newdata dataframe
    #index.newdata <- which(nks.newdata > 0)

    ######################################################
    ## Split data quantities by factor level combinations
    ######################################################

    index.factor.variables <- which(attr(terms(object), "dataClasses") == "factor")
    index.numeric.variables <- which(attr(terms(object), "dataClasses") == "numeric")
    index.data.factor.variables <- which(colnames(data) %in% names(index.factor.variables))
    index.data.numeric.variables <- which(colnames(data) %in% names(index.numeric.variables))    
    factors <- lapply(index.data.factor.variables, function(j) as.numeric(as.factor(data[, j])))

    ## split model matrix by factor level combinations
    split.X.factors <- split(X, factors, drop = FALSE)
    preds.by.factors <- lapply(split.X.factors, function(x){ 
      mat <- matrix(x, ncol = ncol(X))
      colnames(mat) <- colnames(X)
      mat
    })

    ## split model matrix by factor level combinations 
    ## and isolate numeric variables
    numeric.by.factors <- lapply(split.X.factors, function(x){ 
      mat <- matrix(x, ncol = ncol(X))
      colnames(mat) <- colnames(X)
      mat[, index.numeric.variables]
    })

    ## split response by factor level combinations
    split.Y.factors <- split(Y, factors, drop = FALSE)
    resps.by.factors <- lapply(split.Y.factors, as.numeric)

    ## split data frame by factor level combinations
    data.by.factors <- split(data, factors, drop = FALSE)

    ## split bin indices by factor level combinations
    bin.index.by.factors <- lapply(data.by.factors, function(x){ 
      if(nrow(x) == 0) return(0)
      mat <- as.matrix(x[, index.data.numeric.variables[-1]], 
        ncol = length(index.data.numeric.variables[-1]))
      find.index(mat = mat, wn = wn, d = ncol(mat))
    })

    ## get binning structure
    nks <- as.numeric(lengths(resps.by.factors))
 
    ## indices of nonempty factor level combinations in 
    ## data dataframe
    index.data <- which(nks > 0)

    ## create model calls when appropriate
    ## obtain OLS estimate 
    ## calculate important quantities for the gaussian 
    ## distribution
    paraconformal <- NULL
    if(parametric){ 
      sd.res <- NULL
      if(family == "gaussian"){
        m1 <- lm(formula, data = data, x = TRUE)
        betaOLS <- betaMLE <- coefficients(m1)
        sd.res <- summary(m1)$sigma
      }
      ## Get MLEs and plugin interval for gamma distribution 
      betaMLE <- shapeMLE <- rateMLE <- 0
      if(family == "Gamma"){
        betaMLE <- coefficients(object)
        shapeMLE <- as.numeric(gamma.shape(object)[1])
        if(link == "identity"){
          rateMLE <- 1 / (cbind(1, X) %*% betaMLE) * shapeMLE
        }
        if(link == "inverse"){
          rateMLE <- (cbind(1, X) %*% betaMLE) * shapeMLE
        }
        if(link == "log"){
          rateMLE <- (1 / exp(cbind(1, X) %*% betaMLE)) * shapeMLE
        }
      }

      ## Get MLEs and plugin interval for Inverse Gaussian Distribution 
      if(family == "inverse.gaussian"){
        betaMLE <- coefficients(object)
      }

      # upfront quantities to improve speed
      shapeMLE.y <- shapeMLE; rateMLE.y <- rbind(rateMLE, 1); betaMLE.y <- betaMLE
      sd.y <- sd.res
      object.y <- object
      data.y <- matrix(0, nrow = nrow(data) + 1, ncol = ncol(data))
      colnames(data.y) <- colnames(data)
      data.y <- as.data.frame(data.y)

      ## conformal script
      out <- mclapply(1:n.pred, mc.cores = cores, FUN = function(j){

        # factor index of current index
        internal.index <- index.newdata[j]

        # important quantities
        x.variables <- newdata[j, ]
        x <- X.newdata[j, ]
        index.bin.data <- bin.index.by.factors[[internal.index]]
        index.bin <- which(index.bin.data == index.bin.newdata[j])
        Xk <- preds.by.factors[[internal.index]][index.bin, ]
        Yk <- resps.by.factors[[internal.index]][index.bin]
        nk <- length(Yk)

        ## conformal scores (density estimator)
        phatxy <- function(z){

          out <- rateMLE.y <- NULL
          data.y[, colnames(data) %in% respname] <- c(Y, z)
          data.y[, !(colnames(data) %in% respname)] <- rbind(newdata, x.variables)

          if(family == "Gamma"){
            if(link == "identity"){
              m1.y <- glm(formula, data = data.y, family = Gamma(link = identity), 
                control = list(maxit = 1e4))
              shapeMLE.y <- as.numeric(gamma.shape(m1.y)[1])
              rateMLE.y <- 1 / (cbind(1, rbind(Xk, x)) %*% 
                coefficients(m1.y)) * shapeMLE.y
            }
            if(link == "inverse"){
              m1.y <- glm(formula, data = data.y, family = Gamma(link = inverse), 
                control = list(maxit = 1e4))
              shapeMLE.y <- as.numeric(gamma.shape(m1.y)[1])              
              rateMLE.y <- (cbind(1, rbind(Xk, x)) %*% 
                coefficients(m1.y)) * shapeMLE.y
            }
            if(link == "log"){
              m1.y <- glm(formula, data = data.y, family = Gamma(link = log), 
                control = list(maxit = 1e4))
              shapeMLE.y <- as.numeric(gamma.shape(m1.y)[1])              
              rateMLE.y <- (1 / exp(cbind(1, rbind(Xk, x)) %*% 
                coefficients(m1.y))) * shapeMLE.y
            }

            out <- dgamma(c(Yk, z), rate = rateMLE.y, shape = shapeMLE.y)
          }

          if(family == "gaussian"){
            m1.y <- lm(formula, data = data.y[-nrow(data.y), ])
            out <- dnorm(c(Yk, z), mean = as.numeric(cbind(1, rbind(Xk, x)) %*% 
              coefficients(m1.y)), sd = summary(m1.y)$sigma)
          }

          if(family == "inverse.gaussian"){
            m1.y <- glm(formula, data = data.y, family = family, 
              control = list(maxit = 1e4))
            out <- dinvgauss(c(Yk, z), mean = 1 / sqrt(cbind(1, rbind(Xk, x)) %*% 
              coefficients(m1.y)))
          }

          out
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
          suppressWarnings(try.lwr <- try(phatxy(y.lwr), silent = TRUE))
          while(class(try.lwr) == "try-error"){
            y.lwr <- y.lwr + steps * prec
            suppressWarnings(try.lwr <- try(phatxy(y.lwr), silent = TRUE))
            steps <- steps + 1
          }
          steps <- 1          
          flag <- FALSE
          while(rank(phatxy(y.lwr))[nk + 1] >= nk.tilde & flag == FALSE){
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
            while(rank(phatxy(y.lwr))[nk + 1] < nk.tilde & y.lwr < max(Yk)){
              y.lwr <- y.lwr + steps * prec
              steps <- steps + 1
            }
          }

          ## upper line search
          if(y.lwr >= y.upr) y.upr <- 2^sign(max(Yk)) * max(Yk)
          prec <- max( min(diff(sort(Yk[Yk >= y.max]))), precision)
          steps <- 1
          suppressWarnings(try.upr <- try(phatxy(y.upr), silent = TRUE))
          while(class(try.upr) == "try-error"){
            y.upr <- y.upr - steps * prec
            suppressWarnings(try.upr <- try(phatxy(y.upr), silent = TRUE))
            steps <- steps + 1
          }
          steps <- 1
          while(rank(phatxy(y.upr))[nk + 1] >= nk.tilde){
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
          while(rank(phatxy(y.upr))[nk + 1] < nk.tilde){
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

      paraconformal <- do.call(rbind, out)
      colnames(paraconformal) <- c("lwr", "upr")
    }


    # ---- The nonparametric conformal implementation -------
    if(class(h) == "NULL") h <- wn
    nonparaconformal <- NULL
    if(nonparametric){  

      ## the nonparametric conformal prediciton region
      index.pred <- unique(index.newdata)
      nonparaconformal <- mclapply(sort(unique(index.pred)), 
        mc.cores = cores, FUN = function(j){
  
        ## inner loop to get nonparametric conformal prediction region 
        ## for distinct bins within a particular factor level combination
        distinct.bins.within.factors <- lapply(newdata.bin.index.by.factors, unique)[[j]]
        region.by.bin <- lapply(distinct.bins.within.factors, function(k){

          index.bin.data <- bin.index.by.factors[[j]]
          index.bin <- which(index.bin.data == k)
          Xk <- preds.by.factors[[j]][index.bin, ]
          Yk <- resps.by.factors[[j]][index.bin]
          nk <- length(Yk)

          ## initial check for enough data within bin
          endpts <- NULL
          nk.tilde <- floor(alpha * (nk + 1))
          if(nk.tilde <= 1){ 
            endpts <- c(min(Yk), max(Yk))  
          }
          
          ## set up a lower (upper lower bound) and upper bound 
          ## (lower upper bound) to start two line searchs in order to 
          ## construct the parametric conformal prediction region
          if(nk.tilde > 1){

            ## nonparametric density
            phatxy <- function(y){                
              int <- which(unlist(lapply(1:nk, FUN = function(l){
                Yknotl <- Yk[-l]
                Ykl <- Yk[l]
                round(sum(dnorm(Yknotl, mean = y, sd = h) 
                  - dnorm(Yknotl, mean = Ykl, sd = h)), 6)
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
            if(prec == 0) prec <- precision / 2
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
            if(prec == 0) prec <- precision / 2
            while(phatxy(y.upr)){
              y.upr <- y.upr + steps * prec
              steps <- steps + 1
            }

            y.seq <- seq(y.lwr, y.upr, by = precision)
            foo <- unlist(lapply(y.seq, phatxy))
            y.seq <- y.seq[foo]
            breaks <- which(round(diff(y.seq), 
              ceiling(log10(1/precision))) != precision)
            endpts <- c(min(y.seq), max(y.seq))
            if(length(breaks) == 1) endpts <- 
              c(min(y.seq), y.seq[breaks], y.seq[breaks+1],  max(y.seq))
            if(length(breaks) > 1){
              for(l in 1:length(breaks)){
                if(l == 1) endpts <- 
                  c(min(y.seq), y.seq[breaks[l]], y.seq[breaks[l]+1])
                if(l != 1 & l != length(breaks)) endpts <- 
                  c(endpts, y.seq[breaks[l]], y.seq[breaks[l]+1])
                if(l == length(breaks)) endpts <- 
                  c(endpts, y.seq[breaks[l]], y.seq[breaks[l]+1], max(y.seq))
              }
            } 
          }
          endpts
        })
      })
    }

    out = list(paraconformal = paraconformal, 
      nonparaconformal = nonparaconformal)

  }

  return(out)
}

