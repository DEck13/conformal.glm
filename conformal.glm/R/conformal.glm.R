

conformal.glm <- function(object, ..., newdata = NULL, alpha = 0.10, 
  cores = 1, bins = 1, parametric = TRUE, nonparametric = FALSE, 
  method = c("transform", "bin", "both"), intercept = TRUE, h = NULL, 
  precision = 0.005){

  # declare output quantities
  out <- NULL
  transformconf <- NULL 
  paraconfbin <- NULL
  nonparaconfbin <- NULL

  ## some important quantities
  object <- update(object, x = TRUE, y = TRUE)
  call <- object$call
  formula <- call$formula
  fam <- object$family
  family <- fam$family
  stopifnot("glm" %in% class(object))
  stopifnot(family %in% c("Gamma", "gaussian", "inverse.gaussian"))
  link <- fam$link
  data <- object$data

  # ---- Parametric transform conformal -------
  if(method %in% c("transform", "both")){

    ## check for factors. If no factors then call regions function
    if(!any(attr(object$terms, "dataClasses")[-1] == "factor")){

      ## some important quantities
      X <- as.matrix(object$x[, -1])
      colnames(X) <- colnames(object$x)[-1]
      if(is.null(newdata)){ 
        newdata <- data
        respname <- all.vars(formula)[1]
        newdata <- newdata[, !(colnames(data) %in% respname)]
      }
      newdata <- as.matrix(newdata)
      colnames(newdata) <- all.vars(formula)[-1]

      ## model matrix for newdata
      tt <- terms(object)
      Terms <- delete.response(tt)
      X.newdata <- model.matrix(~ ., model.frame(Terms, data.frame(newdata)))[, -1]
      X.newdata <- as.matrix(X.newdata)

      sd.res <- NULL
      betaMLE <- shapeMLE <- rateMLE <- NULL
      ## Get MLEs for Gaussian Distribution 
      if(family == "gaussian"){
        m1 <- lm(formula, data = data, x = TRUE)
        betaMLE <- coefficients(m1)
        sd.res <- summary(m1)$sigma
      }
      ## Get MLEs for gamma distribution 
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
      ## Get MLEs for Inverse Gaussian Distribution 
      if(family == "inverse.gaussian"){
        betaMLE <- coefficients(object)
      }


      ## transform conformal with no factors 
      XnewX <- matrix(0, nrow = n+1, ncol = ncol(X))
      transformconformal <- function(newX){
  
        newX <- matrix(newX, nrow = 1)
        lwr <- upr <- NULL
        rateMLEnewX <- NULL
        ## Get new MLEs and candidate boundary points of 
        ## transform conformal prediction region
        if(family == "gaussian"){
          muX <- (cbind(1, newX) %*% betaMLE)
          lwr <- qnorm(1.25 * alpha, mean = muX, sd = sd.res)
          upr <- qnorm(1 - 1.25 * alpha, mean = muX, sd = sd.res)
        } 
        if(family == "Gamma"){
          if(link == "identity"){
            rateMLEnewX <- 1 / (cbind(1, newX) %*% betaMLE) * shapeMLE
          }
          if(link == "inverse"){
            rateMLEnewX <- (cbind(1, newX) %*% betaMLE) * shapeMLE
          }
          if(link == "log"){
            rateMLEnewX <- (1 / exp(cbind(1, newX) %*% betaMLE)) * shapeMLE
          }
          lwr <- qgamma(1.25 * alpha, rate = rateMLEnewX, shape = shapeMLE)
          upr <- qgamma(1 - 1.25 * alpha, rate = rateMLEnewX, shape = shapeMLE)
        }
        if(family == "inverse.gaussian"){
          betaMLE <- coefficients(object)
        }

        ## augmented data with small upper boundary 
        ## used to get upper boundary of transform conformal prediction region
        aug.data <- rbind(data, 
          c(upr, newX[, colnames(X) %in% all.vars(formula)[-1]])) 
        objectnewX = glm(formula, family = fam, 
          data = aug.data, x = TRUE)
        XnewX <- objectnewX$x
        betaMLEnewX <- coefficients(objectnewX)

        ## transform conformal for linear regression model
        sd.resnewX <- sd.res
        if(family == "gaussian"){
          muX <- XnewX %*% betaMLEnewX
          m1 <- lm(formula, data = aug.data, x = TRUE)
          sd.resnewX <- summary(m1)$sigma
          critval <- c(alpha/2, 1 - alpha/2)

          ## get upper boundary of transform conformal prediction region
          U.upr <- pnorm(upr, mean = muX[n+1, ], sd = sd.resnewX)
          condition <- FALSE
          while(!condition){
            aug.data <- rbind(data, 
              c(upr, newX[, colnames(X) %in% all.vars(formula)[-1]])) 
            objectnewX = lm(formula, data = aug.data, x = TRUE)
            betaMLEnewX <- coefficients(objectnewX)
            XnewX <- objectnewX$x
            munewX <- XnewX %*% betaMLEnewX
            sd.resnewX <- summary(objectnewX)$sigma
            #HDCI <- hdi(qnorm, 1 - alpha, mean = munewX[n+1, 1], sd = sd.resnewX)
            #critval <- pnorm(HDCI, mean = munewX[n+1, 1], sd = sd.resnewX)
            Unorm <- pnorm(aug.data[, 1], mean = munewX, sd = sd.resnewX)
            U.upr <- Unorm[n+1]
            condition <- U.upr - quantile(Unorm, probs = 1 - alpha/2) > 0
            upr <- upr + precision
          }

          ## get lower boundary of transform conformal prediction region
          U.lwr <- pnorm(lwr, mean = muX[n+1, ], sd = sd.res)
          condition <- FALSE
          while(!condition){
            aug.data <- rbind(data, 
              c(lwr, newX[, colnames(X) %in% all.vars(formula)[-1]])) 
            objectnewX = lm(formula, data = aug.data, x = TRUE)
            betaMLEnewX <- coefficients(objectnewX)
            XnewX <- objectnewX$x
            munewX <- XnewX %*% betaMLEnewX
            sd.resnewX <- summary(objectnewX)$sigma
            #HDCI <- hdi(qnorm, 1 - alpha, mean = munewX[n+1, 1], sd = sd.resnewX)
            #critval <- pnorm(HDCI, mean = munewX[n+1, 1], sd = sd.resnewX)
            Unorm <- pnorm(aug.data[, 1], mean = munewX, sd = sd.resnewX)
            U.lwr <- Unorm[n+1]
            condition <- quantile(Unorm, probs = alpha/2) - U.lwr > 0
            lwr <- lwr - precision
          }
        }

        ## transform conformal for gamma regression model
        if(family == "Gamma"){
          shapeMLEnewX <- as.numeric(gamma.shape(objectnewX)[1])
          if(link == "identity"){
            rateMLEnewX <- 1 / (XnewX %*% betaMLEnewX) * shapeMLEnewX
          }
          if(link == "inverse"){
            rateMLEnewX <- (XnewX %*% betaMLEnewX) * shapeMLEnewX
          }
          if(link == "log"){
            rateMLEnewX <- (1 / exp(XnewX %*% betaMLEnewX)) * shapeMLEnewX
          }        

          ## get percentiles for highest density region
          HDCI <- do.call(rbind, 
            lapply(1:nrow(aug.data), function(j){ 
            hdi(qgamma, 1 - alpha, shape = shapeMLEnewX, rate = rateMLEnewX[j, 1])
          }))
          critval <- pgamma(HDCI[n+1, ], shape = shapeMLEnewX, rate = rateMLEnewX[n+1, 1])

          ## get upper boundary of transform conformal prediction region
          U.upr <- pgamma(upr, rate = rateMLEnewX[n+1, ], shape = shapeMLEnewX)
          condition <- FALSE
          while(!condition){
            aug.data <- rbind(data, 
              c(upr, newX[, colnames(X) %in% all.vars(formula)[-1]])) 
            objectnewX = glm(formula, family = fam, data = aug.data, x = TRUE)
            betaMLEnewX <- coefficients(objectnewX)
            XnewX <- objectnewX$x
            shapeMLEnewX <- as.numeric(gamma.shape(objectnewX)[1])
            if(link == "identity"){
              rateMLEnewX <- 1 / (XnewX %*% betaMLEnewX) * shapeMLEnewX
            }
            if(link == "inverse"){
              rateMLEnewX <- (XnewX %*% betaMLEnewX) * shapeMLEnewX
            }
            if(link == "log"){
              rateMLEnewX <- (1 / exp(XnewX %*% betaMLEnewX)) * shapeMLEnewX
            }
            HDCI <- hdi(qgamma, 1 - alpha, shape = shapeMLEnewX, rate = rateMLEnewX[n+1, 1])
            critval <- pgamma(HDCI, shape = shapeMLEnewX, rate = rateMLEnewX[n+1, 1])
            Ugamma <- pgamma(aug.data[, 1], rate = rateMLEnewX, shape = shapeMLEnewX)
            U.upr <- Ugamma[n+1]
            condition <- U.upr - quantile(Ugamma, probs = critval[2]) > 0
            upr <- upr + precision
          }

          ## get lower boundary of transform conformal prediction region
          U.lwr <- pgamma(lwr, rate = rateMLEnewX[n+1, ], shape = shapeMLEnewX)
          condition <- FALSE
          while(!condition){
            aug.data <- rbind(data, 
              c(lwr, newX[, colnames(X) %in% all.vars(formula)[-1]])) 
            objectnewX = glm(formula, family = fam, data = aug.data, x = TRUE)
            betaMLEnewX <- coefficients(objectnewX)
            XnewX <- objectnewX$x
            shapeMLEnewX <- as.numeric(gamma.shape(objectnewX)[1])
            if(link == "identity"){
              rateMLEnewX <- 1 / (XnewX %*% betaMLEnewX) * shapeMLEnewX
            }
            if(link == "inverse"){
              rateMLEnewX <- (XnewX %*% betaMLEnewX) * shapeMLEnewX
            }
            if(link == "log"){
              rateMLEnewX <- (1 / exp(XnewX %*% betaMLEnewX)) * shapeMLEnewX
            }
            HDCI <- hdi(qgamma, 1 - alpha, shape = shapeMLEnewX, rate = rateMLEnewX[n+1, 1])
            critval <- pgamma(HDCI, shape = shapeMLEnewX, rate = rateMLEnewX[n+1, 1])
            Ugamma <- pgamma(aug.data[, 1], rate = rateMLEnewX, shape = shapeMLEnewX)
            U.lwr <- Ugamma[n+1]
            condition <- quantile(Ugamma, probs = critval[1]) - U.lwr > 0
            lwr <- lwr - precision
            if(lwr < 0){
              lwr <- 0.001
              break
            }
          }
        }
        c(lwr, upr)
      }

      transformconf <- do.call(rbind, mclapply(1:nrow(X.newdata), mc.cores = cores, 
        FUN = function(j){ transformconformal(X.newdata[j, ]) }))

    }

    if(any(attr(object$terms, "dataClasses")[-1] == "factor")){

      ## some important quantities    
      X <- object$x[, -1]
      Y <- object$y
      n <- length(Y)
      p <- ncol(X) - 1
      respname <- all.vars(formula)[1]
      terms <- c(respname, attr(terms(object), "term.labels"))
      data <- data[, colnames(data) %in% terms]
      d <- ncol(data)
      wn <- min(1/ floor(1 / (log(n)/n)^(1/(d+3))), 1/2)
      if(class(bins) != "NULL") wn <- 1 / bins

      #newdata <- NULL
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

      #####################################
      ## transform conformal with factors #
      #####################################

    }



    out = list(transformconf = transformconf, 
      paraconfbin = paraconfbin, 
      nonparaconfbin = nonparaconfbin)
  }


  

  # ---- Conformal prediction w/ binning -------
  if(method %in% c("bin", "both")){
  
    ## check for factors. If no factors then call regions function
    if(!any(attr(object$terms, "dataClasses")[-1] == "factor")){

      ## some important quantities
      if(is.null(newdata)){ 
        newdata <- data
        respname <- all.vars(formula)[1]
        newdata <- newdata[, !(colnames(data) %in% respname)]
      }
      newdata <- as.matrix(newdata)


      int <- regions(formula = formula, data = data, 
        newdata = newdata, family = family, link = link, 
        alpha = alpha, cores = cores, bins = bins, 
        parametric = parametric, nonparametric = nonparametric, 
        h = h, precision = precision)

      paraconfbin <- int$paraconformal
      nonparaconfbin <- int$nonparaconformal
      out = list(transformconf = transformconf, 
        paraconfbin = paraconfbin, 
        nonparaconfbin = nonparaconfbin)    
    }

    if(any(attr(object$terms, "dataClasses")[-1] == "factor")){

      ## some important quantities
      X <- object$x[, -1]
      Y <- object$y
      n <- length(Y)
      p <- ncol(X) - 1
      respname <- all.vars(formula)[1]
      terms <- c(respname, attr(terms(object), "term.labels"))
      data <- data[, colnames(data) %in% terms]
      d <- ncol(data)
      wn <- min(1/ floor(1 / (log(n)/n)^(1/(d+3))), 1/2)
      if(class(bins) != "NULL") wn <- 1 / bins

      #newdata <- NULL
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
        mat <- as.matrix(x[, index.newdata.numeric.variables], 
          ncol = length(index.newdata.numeric.variables))
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

      # ---- The parametric conformal implementation w/ binning -------
      ## create model calls when appropriate
      ## obtain OLS estimate 
      ## calculate important quantities for the gaussian 
      ## distribution
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
        shapeMLE.y <- shapeMLE 
        rateMLE.y <- rbind(rateMLE, 1) 
        betaMLE.y <- betaMLE
        sd.y <- sd.res
        object.y <- object
        #data.y <- matrix(0, nrow = nrow(data) + 1, ncol = ncol(data))
        #colnames(data.y) <- colnames(data)
        #data.y <- as.data.frame(data.y)
        #data.y <- data
        #data.y <- rbind(data.y, 0)
        data.y <- rbind(data, data[1, ])

        ## conformal script
        out <- mclapply(1:n.pred, mc.cores = cores, FUN = function(j){

          # factor index of current index
          internal.index <- index.newdata[j]

          # important quantities
          xnew <- newdata[j, ]
          xnew.modmat <- X.newdata[j, ]
          index.bin.data <- bin.index.by.factors[[internal.index]]
          index.bin <- which(index.bin.data == index.bin.newdata[j])
          Xk <- preds.by.factors[[internal.index]][index.bin, ]
          Yk <- resps.by.factors[[internal.index]][index.bin]
          nk <- length(Yk)

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
            suppressWarnings(try.lwr <- try(density_scores(y.lwr), silent = TRUE))
            while(class(try.lwr) == "try-error"){
              y.lwr <- y.lwr + steps * prec
              suppressWarnings(try.lwr <- try(phatxy(y.lwr), silent = TRUE))
              steps <- steps + 1
            }
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
            suppressWarnings(try.upr <- try(density_scores(y.upr), silent = TRUE))
            while(class(try.upr) == "try-error"){
              y.upr <- y.upr - steps * prec
              suppressWarnings(try.upr <- try(density_scores(y.upr), silent = TRUE))
              steps <- steps + 1
            }
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

        paraconfbin <- do.call(rbind, out)
        colnames(paraconfbin) <- c("lwr", "upr")
      }

      # ---- The nonparametric conformal implementation -------
      if(class(h) == "NULL") h <- wn
      if(nonparametric){  

        ## the nonparametric conformal prediciton region
        index.pred <- unique(index.newdata)
        nonparaconfbin <- mclapply(sort(unique(index.pred)), 
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
              kernxy <- function(y){                
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
              while(!kernxy(y.lwr)){
                y.lwr <- y.lwr + steps * prec
                steps <- steps + 1
              }
              steps <- 1
              prec <- min( min(diff(sort(Yk[Yk <= quant.Yk[1]]))), precision)
              if(prec == 0) prec <- precision / 2
              while(kernxy(y.lwr)){
                y.lwr <- y.lwr - steps * prec
                steps <- steps + 1    
              }

              # upper line search
              y.upr <- y.max <- max(Yk)
              prec <- max( min(diff(sort(Yk[Yk >= quant.Yk[2]]))), precision)
              steps <- 1
              while(!kernxy(y.upr)){
                y.upr <- y.upr - steps * prec
                steps <- steps + 1        
              }
              steps <- 1
              prec <- min( min(diff(sort(Yk[Yk >= quant.Yk[2]]))), precision)
              if(prec == 0) prec <- precision / 2
              while(kernxy(y.upr)){
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
    }    
    out = list(transformconf = transformconf, 
      paraconfbin = paraconfbin, 
      nonparaconfbin = nonparaconfbin)
  }
  return(out)
}








