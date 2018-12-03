
find.index <- function(mat, wn, k){
  A <- seq(from = 0, to = 1 - wn, by = wn)
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
  alpha = 0.10, cores = 1, bins = NULL, intercept = TRUE, parametric = TRUE, 
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
  key <- cbind(X, Y, index)

  #subkey <- lapply(unique(index.pred), FUN = function(j){
  #  datak <- key[key[, p + 2] == j, ]
  #  datak
  #})
  #subkey <- mclapply(sort(unique(index)), FUN = function(j){
  #  datak <- key[key[, p + 2] == j, ]
  #  datak
  #}, mc.cores = cores)
  #subkey <- split(cbind(X, Y), f = as.factor(index))

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
        #index.j <- which(index.pred[j] == indices.pred)
        #datak <- matrix(subkey[[index.j]], ncol = p + 2)
        #Xk <- matrix(datak[, 1:p], ncol = p)
        #Yk <- datak[, p+1]
        Yk <- key[key[, p + 2] == index.pred[j], ][, p+1]
        nk <- length(Yk)

        ## conformal scores
        phatxy <- function(z){
          out <- NULL
          data.y[, colnames(data) %in% respname] <- c(Y, z)
          data.y[, !(colnames(data) %in% respname)] <- rbind(X.variables, x.variables)
          data.y <- as.data.frame(data.y)

          if(family == "Gamma"){
            m1.y <- glm(formula, data = data.y, family = family)
            shapeMLE.y <- as.numeric(gamma.shape(m1.y)[1])
            if(link == "identity"){
              rateMLE.y <- 1 / (cbind(1, x) %*% coefficients(m1.y)) * shapeMLE.y
            }
            if(link == "inverse"){
              rateMLE.y <- cbind(1, x) %*% coefficients(m1.y) * shapeMLE.y
            }
            if(link == "log"){
              rateMLE.y <- (1 / exp(cbind(1, x) %*% coefficients(m1.y))) * shapeMLE.y
            }

            out <- dgamma(c(Y, z), 
              rate = cbind(1, rbind(X, x)) %*% coefficients(m1.y) * shapeMLE.y, 
              shape = shapeMLE.y)
          }

          if(family == "gaussian"){
            m1.y <- lm(formula, data = data.y)
            out <- dnorm(c(Y, z), 
              mean = as.numeric(cbind(1, rbind(X, x)) %*% coefficients(m1.y)), 
              sd = summary(m1.y)$sigma)
          }

          if(family == "inverse.gaussian"){
            m1.y <- glm(formula, data = data.y, family = family)
            out <- dinvgauss(c(Y, z), mean = 1 / sqrt(cbind(1, rbind(X, x)) %*% coefficients(m1.y)))
          }

          out
        }

        ## set up range of candidate points used to 
        ## construct the parametric conformal prediction 
        ## region
        lwr <- min(Yk)
        upr <- max(Yk)
        if(lwr > 0) lwr <- lwr / 1.5
        if(lwr < 0) lwr <- lwr * 1.5
        if(upr > 0) upr <- upr * 1.5
        if(upr < 0) upr <- upr / 1.5

        ## an intial crude approximation of the 
        ## parametric conformal prediction region 
        nk.tilde <- floor(alpha * (nk + 1))
        y.lwr <- lwr
        y.upr <- upr

        if(nk.tilde > 0){

          ## perform a crude crude search 
          prec <- max( min(diff(sort(Yk))), 0.001)
          crudewidth <- sqrt((upr - lwr)*prec*2*alpha)
          crude.seq.y <- seq(from = lwr, to = upr, 
            by = crudewidth)
          quant.Yk <- quantile(Yk, probs = c(alpha, 1 - alpha))
          crude.seq.y <- as.list(
            crude.seq.y[!(crude.seq.y > quant.Yk[1] & crude.seq.y < quant.Yk[2])])
          crude.search <- lapply(crude.seq.y, 
            FUN = function(y){ 
              foo <- cbind(phatxy(y), c(index, index.pred[j]))
              bar <- foo[which(foo[, 2] == index.pred[j]), -2]
              baz <- rank(bar)[nk + 1]
              baz
            })
          
          ## get nk threshold        
          cand <- which(crude.search >= nk.tilde)
          breaks <- diff(cand) ## corresponds to number of modes

          if(unique(breaks) == 1){ ## modes occur in an interval 
            endpt1.lwr <- endpt2.lwr <- endpt1.upr <- 
              endpt2.upr <- 0
            if(min(cand) == 1){ 
              endpt1.lwr <- lwr - crudewidth
              endpt2.lwr <- lwr
              if(family != "gaussian") endpt1.lwr <- max(endpt1.lwr, 0.00001)
            }
            if(max(cand) == length(crude.seq.y)){ 
              endpt2.upr <- upr + crudewidth
              endpt1.upr <- upr
            }
            if(min(cand) > 1){
              endpt1.lwr <- min(crude.seq.y[[min(cand) - 1]], 
                crude.seq.y[[min(cand)]],
                crude.seq.y[[min(cand) + 1]])
              endpt2.lwr <- max(crude.seq.y[[min(cand) - 1]], 
                crude.seq.y[[min(cand)]], 
                crude.seq.y[[min(cand) + 1]])
            }
            if(max(cand) < length(crude.seq.y)){ 
              endpt2.upr <- max(crude.seq.y[[max(cand) + 1]], 
                crude.seq.y[[max(cand)]], 
                crude.seq.y[[max(cand) - 1]])
              endpt1.upr <- min(crude.seq.y[[max(cand) + 1]], 
                crude.seq.y[[max(cand)]],
                crude.seq.y[[max(cand) - 1]])
            }  

            if(family == "Gamma") if(lwr <= 0.01) endpt1.lwr <- 0.00001

            precise.seq.y1 <- as.list(seq(from = endpt1.lwr, 
              to = endpt2.lwr, length = floor(crudewidth/prec)))
            precise.search1 <- lapply(precise.seq.y1, 
              FUN = function(y){
                foo <- cbind(phatxy(y), c(index, index.pred[j]))
                bar <- foo[which(foo[, 2] == index.pred[j]), -2]
                baz <- rank(bar)[length(bar)]
                baz
              })
            cand1 <- min(which(precise.search1 >= nk.tilde))
            y.lwr <- precise.seq.y1[[cand1]]

            precise.seq.y2 <- as.list(seq(from = endpt1.upr, 
              to = endpt2.upr, length = floor(crudewidth/prec)))
            precise.search2 <- lapply(precise.seq.y2, 
              FUN = function(y){ 
                foo <- cbind(phatxy(y), c(index, index.pred[j]))
                bar <- foo[which(foo[, 2] == index.pred[j]), -2]
                baz <- rank(bar)[length(bar)]
                baz
              })
            cand2 <- max(which(precise.search2 >= nk.tilde))
            y.upr <- precise.seq.y2[[cand2]]
          }
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
        x <- newdata[j, ]
        index.j <- which(index.pred[j] == indices.pred)
        datak <- matrix(key[key[, p + 2] == index.pred[j], ][, -c(p+2)], ncol = p+1)
        Xk <- matrix(datak[, 1:p], ncol = p)
        Yk <- datak[, p+1]
        nk <- nrow(datak)

        ## nonparametric density
        phatxy <- function(y){                
          out <- which(unlist(lapply(1:nk, FUN = function(j){
            Yknotj <- Yk[-j]
            Ykj <- Yk[j]
            sum(dnorm(Yknotj, mean = y, sd = hn) 
              - dnorm(Yknotj, mean = Ykj, sd = hn))
          })) >= 0)
          if(length(out) == 0) out <- -1
          out
        }

        ## set up range of candidate points used to 
        ## construct the parametric conformal prediction 
        ## region
        lwr <- min(Yk)
        upr <- max(Yk)
        if(lwr > 0) lwr <- lwr / 1.5
        if(lwr < 0) lwr <- lwr * 1.5
        if(upr > 0) upr <- upr * 1.5
        if(upr < 0) upr <- upr / 1.5

        ## an intial crude approximation of the 
        ## parametric conformal prediction region 
        nk.tilde <- floor(alpha * (nk + 1))
        #crudewidth <- sqrt((upr - lwr)*0.001) # chosen to minimze 
           # the number of candidate points to search over when we 
           # desire a final precision of 0.001
        prec <- max( min(diff(sort(Yk))), 0.001)   
        crudewidth <- sqrt((upr - lwr)*prec*2*alpha)                    
        #crude.seq.y <- seq(from = lwr, to = upr, by = crudewidth)
        #crude.search <- apply(matrix(crude.seq.y), 1, 
        #  FUN = function(y){
        #    int <- out <- phatxy(y)
        #    out <- length(int)
        #    if(out == 1){
        #      if(int == -1) out <- 0
        #    }
        #    out
        #  })
        #crude.seq.y <- as.list(seq(from = lwr, to = upr, by = crudewidth))
        crude.seq.y <- seq(from = lwr, to = upr, 
          by = crudewidth)
        quant.Yk <- quantile(Yk, probs = c(alpha, 1 - alpha))
        crude.seq.y <- as.list(
          crude.seq.y[!(crude.seq.y > quant.Yk[1] & crude.seq.y < quant.Yk[2])])
        crude.search <- lapply(crude.seq.y, FUN = function(y){
            int <- out <- phatxy(y)
            out <- length(int)
            if(out == 1){
              if(int == -1) out <- 0
            }
            out
          })         
        cand <- which(crude.search >= nk.tilde)
        breaks <- diff(cand) ## corresponds to number of modes

        y.lwr <- y.upr <- 0
        if(unique(breaks) == 1){ ## modes occur in an interval 
          endpt1.lwr <- endpt2.lwr <- endpt1.upr <- 
            endpt2.upr <- 0
          if(min(cand) == 1){ 
            endpt1.lwr <- lwr - crudewidth 
            endpt2.lwr <- lwr
          }
          if(max(cand) == length(crude.seq.y)){ 
            endpt2.upr <- upr + crudewidth
            endpt1.upr <- upr
          }
          if(min(cand) > 1){
            endpt1.lwr <- crude.seq.y[[min(cand) - 1]] 
            endpt2.lwr <- crude.seq.y[[min(cand)]] 
          }
          if(max(cand) < length(crude.seq.y)){ 
            endpt2.upr <- crude.seq.y[[max(cand) + 1]]
            endpt1.upr <- crude.seq.y[[max(cand)]]
          }  

          precise.seq.y1 <- as.list(seq(from = endpt1.lwr, 
            to = endpt2.lwr, length = floor(crudewidth/prec)))
          precise.search1 <- lapply(precise.seq.y1, FUN = function(y){
              int <- out <- phatxy(y)
              out <- length(int)
              if(out == 1){
                if(int == -1) out <- 0
              }
              out
            }) 
          cand1 <- min(which(precise.search1 >= nk.tilde))
          y.lwr <- precise.seq.y1[[cand1]]
 
          precise.seq.y2 <- as.list(seq(from = endpt1.upr, 
            to = endpt2.upr, length = floor(crudewidth/prec)))
          precise.search2 <- lapply(precise.seq.y2, FUN = function(y){
              int <- out <- phatxy(y)
              out <- length(int)
              if(out == 1){
                if(int == -1) out <- 0
              }
              out
            }) 
          cand2 <- max(which(precise.search2 >= nk.tilde))
          y.upr <- precise.seq.y2[[cand2]]
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


