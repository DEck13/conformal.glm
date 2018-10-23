

regions <- function(formula, data, newdata, family = "gaussian", link, 
  alpha = 0.10, cores = 6, bins = NULL, intercept = TRUE, parametric = TRUE, 
  LS = FALSE, nonparametric = FALSE){

  ## initial quantities
  respname <- all.vars(formula)[1]
  Y <- data[, colnames(data) %in% respname]
  n <- length(Y)
  n.pred <- nrow(newdata)
  convergence <- 0
 
  ## create model calls when appropriate
  ## obtain OLS estimate 
  ## calculate important quantities for the gaussian 
  ## distribution 
  m1 <- lm(formula, data = data, x = TRUE)
  X <- matrix(m1$x[, -1], nrow = n)
  X.variables <- as.matrix(data[,-which(colnames(data) %in% respname)], 
    nrow = n)
  k <- ncol(X.variables)
  betaOLS <- betaMLE <- m1$coefficients
  p <- length(betaOLS) - 1
  InvFish <- p1 <- pred <- sepred <- #interval.glm <- 
    interval.plugin <- NULL
  sd.res <- sqrt(1/(n-p) * sum(m1$residuals^2))

  ## Get plugin interval for Gaussian Distribution 
  if(family == "gaussian"){
    p1 <- predict(m1, newdata = data.frame(newdata), se.fit = T)
    pred <- p1$fit
    sepred <- p1$se.fit * sqrt(n)

    ## compute glm prediction region
    #interval.plugin <- cbind(pred - qnorm(1 - alpha/2) * sepred, 
    #  pred + qnorm(1 - alpha/2) * sepred)
    #sepred <- sd.res * sqrt(1 + 1/n) ## more needed
    interval.plugin <- cbind(pred - qnorm(1 - alpha/2) * sepred, 
      pred + qnorm(1 - alpha/2) * sepred)
  }


  ## Get MLEs and plugin interval for gamma gamma distribution 
  shapeMLE <- rateMLE <- 0
  if(family == "Gamma"){
    m1 <- glm(formula, data = data, family = "Gamma")
    InvFish <- vcov(m1)
    betaMLE <- m1$coefficients
    shapeMLE <- as.numeric(gamma.shape(m1)[1])
    rateMLE <- cbind(1, X) %*% betaMLE * shapeMLE

    p1 <- predict(m1, type = "response", 
      newdata = data.frame(newdata), se.fit = TRUE)
    pred <- p1$fit
    #sepred <- p1$se.fit
    ## compute glm prediction region 
    #interval.glm <- cbind(pred - qnorm(1 - alpha/2) * sepred, 
    #  pred + qnorm(1 - alpha/2) * sepred)

    ## compute interval prediction region 
    sepred <- sqrt(1 / shapeMLE * pred^2)
    interval.plugin <- cbind(pred - qnorm(1 - alpha / 2) * sepred, 
      pred + qnorm(1 - alpha / 2) * sepred)
  }


  ## Get MLEs and plugin interval for Inverse Gaussian Distribution 
  if(family == "inverse.gaussian"){
    m1 <- glm(formula, data = data, family = family)
    betaMLE <- m1$coefficients
    pred <- predict(m1, type = "response", newdata = data.frame(newdata))
    gx <- 1 / sqrt( cbind(1, X) %*% betaMLE  )
    scaleMLE <- 1 / mean((Y - gx)^2 / (Y * gx^2) )
    sepred <- sqrt(pred^3 * scaleMLE)
    
    ## compute glm prediction region 
    interval.plugin <- cbind(pred - qnorm(1 - alpha/2) * sepred, 
      pred + qnorm(1 - alpha/2) * sepred)
  }

  ## set up partition
  #wn <- min(1 / ceiling((n/log(n))^(k/(k+3))), 
  #  1/3) # width of partition (max width = 1/3)
  #wn <- min(1 / floor((n/log(n))^(1/2)), 
  #  1/3) # width of partition (max width = 1/3)
  #wn <- min(1 / floor(n / 100), 1/3)
  wn <- min(1/ floor(1 / (log(n)/n)^(1/(k+3))), 1/2)
  #if(n <= 200) wn <- 1/2
  if(class(bins) != "NULL") wn <- 1 / bins
  A1 <- seq(from = 0, to = 1 - wn, by = wn)
  A2 <- seq(from = wn, to = 1, wn)
  A <- cbind(A1, A2)

  ## function to allocate X and Y data to all relevant 
  ## partitions via indexing
  ## takes an nXk matrix as an argument where k is the 
  ## number of unique variables
  findindex <- function(mat){
    n.out <- nrow(mat)
    out <- rep(0, n.out)
    for(j in 1:n.out){
      foo <- 0
      for(i in 1:k){
        if(i < k){
          foo <- foo + 
            (max(which(A[, 1] < mat[j, i])) - 1) * (1/wn)^(k-i)
        }
        if(i == k){
          foo <- foo + 
            max(which(A[, 1] < mat[j, i]))
        }
      }
      out[j] <- foo
    }
    out
  }

  ## important internal quantities for our prediction regions   
  ## create seperate data frames with respect to each partition 
  ## ignores the intercept of the newdata matrix
  index <- findindex(X.variables)
  index.pred <- findindex(newdata)

  key <- cbind(X, Y, index)
  #subkey <- lapply(unique(index.pred), FUN = function(j){
  #  datak <- key[key[, p + 2] == j, ]
  #  datak
  #})
  subkey <- mclapply(sort(unique(index)), FUN = function(j){
    datak <- key[key[, p + 2] == j, ]
    datak
  }, mc.cores = cores)
  #subkey <- split(cbind(X, Y), f = as.factor(index))

  ## get newdata in form that is accepted by LS conformal 
  ## and our paramteric conformal implementation
  #newform <- update(formula, NULL ~ ., data = data)
  newdata.variables <- as.matrix(
    model.matrix(~ ., data.frame(newdata))[, -1])
  indices.pred <- sort(unique(index.pred))

  paraconformal <- LSconformal <- nonparaconformal <- NULL
  if(parametric == TRUE){
    # upfront quantities to improve speed
    #newdata.list <- lapply(1:nrow(newdata), 
    #  FUN = function(j) newdata[j, ])

    # ---- The parametric conformal implementation -------
    ## Parametric conformal prediction region for each 
    ## desired predictor combination 
    ## ignores the intercept of the newdata matrix
    Copt <- function(newdata){ 
      out <- mclapply(1:n.pred, mc.cores = cores, FUN = function(j){
        #x <- newdata.list[[j]]
        x <- newdata[j, ]
        index.j <- which(index.pred[j] == indices.pred)
        #datak <- subkey[[index.j]], ncol = k+1)
        datak <- matrix(subkey[[index.j]], ncol = p + 2)
        Xk <- matrix(datak[, 1:p], ncol = p)
        Yk <- datak[, p+1]
        nk <- nrow(datak)
        rateMLE <- cbind(1, Xk) %*% betaMLE * shapeMLE
        rate.x <- as.numeric(c(1, x) %*% betaMLE * shapeMLE)
        lam1 <- nk / (nk + 1)
        lam2 <- 1 - lam1
        phatxy <- function(y){
          out <- NULL

          if(family == "Gamma"){
            ## phat|Ak
            out <- lam1 * mean(dgamma(y, shape = shapeMLE, 
              rate = rateMLE)) + 
            ## phat|x
              lam2 * dgamma(y, shape = shapeMLE, rate = rate.x)  
          }

          if(family == "gaussian"){
            ## phat|Ak
            out <- lam1 * mean(dnorm(y, 
              mean = cbind(1, Xk) %*% betaMLE, sd = sd.res)) + 
            ## phat|x
            lam2 * dnorm(y, mean = c(1, x) %*% betaMLE, sd = sd.res)    
          }

          if(family == "inverse.gaussian"){
            out <- lam1 * mean(dinvgauss(y, 
              mean = 1 / sqrt(cbind(1, Xk) %*% betaMLE))) + 
            ## phat|x
            lam2 * dinvgauss(y, 
              mean = 1 / sqrt(c(1, x) %*% betaMLE))
          }

          out
        }



        ## set up range of candidate points used to 
        ## construct the parametric conformal prediction 
        ## region
        lwr <- min(Yk)
        if(min(Yk) < 0) lwr <- 2 * lwr
        if(min(Yk) >= 0) lwr <- (1/2) * lwr
        upr <- max(Yk)
        if(max(Yk) < 0) upr <- (1/2) * upr
        if(max(Yk) >= 0) upr <- 2 * upr

        ## an intial crude approximation of the 
        ## parametric conformal prediction region 
        nk.tilde <- floor(alpha * (nk + 1))
        y.lwr <- lwr
        y.upr <- upr

        if(nk.tilde > 0){

          ## calculate phatxy for all observed responses
          phatxyY <- lapply(as.list(Yk), FUN = phatxy)
          #phatxyY <- phatxy(Yk)

          ## perform a crude crude search 
          crudewidth <- sqrt((upr - lwr)*0.001)
          crude.seq.y <- as.list(seq(from = lwr, to = upr, 
            by = crudewidth))
          crude.search <- lapply(crude.seq.y, 
            FUN = function(y){ 
              sum(phatxy(y) >= phatxyY)
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

            precise.seq.y1 <- as.list(seq(from = endpt1.lwr, 
              to = endpt2.lwr, length = floor(crudewidth/0.001)))
            precise.search1 <- lapply(precise.seq.y1, 
              FUN = function(y) sum(phatxy(y) >= phatxyY))
            #precise.search1 <- apply(
            #  matrix(precise.seq.y1), 1, FUN = function(y){
            #    length(which(phatxy(y) >= phatxyY))
            #}) 
            cand1 <- min(which(precise.search1 >= nk.tilde))
            y.lwr <- precise.seq.y1[[cand1]]

            #precise.seq.y2 <- seq(from = endpt1.upr, 
            #  to = endpt2.upr, by = 0.001)
            #precise.search2 <- apply(matrix(precise.seq.y2), 
            #  1, FUN = function(y){
            #    length(which(phatxy(y) >= phatxyY))
            #  })
            precise.seq.y2 <- as.list(seq(from = endpt1.upr, 
              to = endpt2.upr, length = floor(crudewidth/0.001)))
            precise.search2 <- lapply(precise.seq.y2, 
              FUN = function(y) sum(phatxy(y) >= phatxyY))
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
  if(LS == TRUE){

    ## Get regression based conformal prediction region 
    ## using conformal.pred in the conformalInference 
    ## package 
    funs <- lm.funs(intercept = intercept)
    train.fun <- funs$train.fun
    predict.fun <- funs$predict.fun
    p1.tibs <- conformal.pred(x = X, y = Y, x0 = newdata.variables, 
      train.fun = train.fun, predict.fun = predict.fun, 
      alpha = alpha, grid.method ="linear",
      num.grid.pts = 999)
    LSconformal <- cbind(p1.tibs$lo, p1.tibs$up)

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
        #datak <- subkey[[index.j]], ncol = k+1)
        datak <- matrix(subkey[[index.j]], ncol = p + 2)
        Xk <- matrix(datak[, 1:p], ncol = p)
        Yk <- datak[, p+1]
        nk <- nrow(datak)
        #x <- newdata[j, ]
        #index.j <- index.pred[j]
        #datak <- subkey[[index.j]]
        #Xk <- datak[, 1:p]
        #Yk <- datak[, p+1]
        #nk <- nrow(datak)

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
        if(min(Yk) < 0) lwr <- 2 * lwr
        if(min(Yk) >= 0) lwr <- (1/2) * lwr
        upr <- max(Yk)
        if(max(Yk) < 0) upr <- (1/2) * upr
        if(max(Yk) >= 0) upr <- 2 * upr

        ## an intial crude approximation of the 
        ## parametric conformal prediction region 
        nk.tilde <- floor(alpha * (nk + 1))
        crudewidth <- sqrt((upr - lwr)*0.001) # chosen to minimze 
           # the number of candidate points to search over when we 
           # desire a final precision of 0.001                    
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
        crude.seq.y <- as.list(seq(from = lwr, to = upr, by = crudewidth))
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
            to = endpt2.lwr, length = floor(crudewidth/0.001)))
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
            to = endpt2.upr, length = floor(crudewidth/0.001)))
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

  #h <- ((log(n)/n)^(1/(1*(p+2)+1)))
  #COPS <- function(x){

  #  index <- index.pred[x]
  #  datak <- subkey[[index]]
  #  Xk <- datak[, 1:p]
  #  Yk <- datak[, p+1]
  #  nk <- nrow(datak)

  #  alpha <- floor((nk+1)*alpha) / (nk+1)

  #  dens <- function(v){
  #    mean(dnorm(Yk, mean = v, sd = h))
  #  }

  #  ## hat(p)^(x,y)(v|Ak)
  #  c1 <- nk / (nk + 1)
  #  dens.new <- function(v, y){
  #    c1 * dens(v) +  (1 - c1) * dnorm(y, mean = v, sd = h)
  #  }

  #  pi <- function(y){
  #    mean(apply(matrix(c(Yk, y)), 1, 
  #      FUN = function(x){ 
  #        as.numeric( dens.new(x, y) <= dens.new(y, y) )
  #    }))
  #  }

  #  a <- max( abs(min(Yk)), abs(max(Yk))  )
  #  y.course <- seq(from = -1.5*a, to = 1.5*a, length = 15)
  #  b <- apply(matrix(y.course), 1, pi)
  #  c <- which(b >= alpha)
  #  y.fine <- seq(from = y.course[min(c)-1], 
  #    to = y.course[min(c)], length = 250)
  #  b2 <- apply(matrix(y.fine), 1, pi)
  #  y.fine2 <- seq(from = y.course[max(c)], 
  #    to = y.course[max(c)+1], length = 250)
  #  b3 <- apply(matrix(y.fine2), 1, pi)

  #  c(y.fine[min(which(b2 >= alpha))], 
  #    y.fine2[max(which(b3 >= alpha))])

  #}

  ### obtain a parametric conformal prediction region for each 
  ### element of the newdata dataframe
  #nonparaconformal <- matrix(0, nrow = nrow(newdata), ncol = 2)
  #for(j in 1:nrow(newdata)){
  #  nonparaconformal[j, ] <- COPS(j)
  #}


    nonparaconformal <- COPS(newdata)

  }

  out = list(paraconformal = paraconformal, 
    nonparaconformal = nonparaconformal,
    LSconformal = LSconformal, #interval.glm = interval.glm, 
    interval.plugin = interval.plugin)
  return(out)
}


