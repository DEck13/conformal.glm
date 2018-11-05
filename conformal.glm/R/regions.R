

regions <- function(formula, data, newdata, family = "gaussian", link, 
  alpha = 0.10, cores = 1, bins = NULL, intercept = TRUE, parametric = TRUE, 
  nonparametric = FALSE){

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
  #InvFish <- p1 <- pred <- sepred <- #interval.glm <- 
  #  interval.plugin <- NULL
  sd.res <- summary(m1)$sigma

  ## Get plugin interval for Gaussian Distribution 
  #if(family == "gaussian"){
  #  p1 <- predict(m1, newdata = data.frame(newdata), se.fit = T)
  #  pred <- p1$fit
  #  sepred <- p1$se.fit * sqrt(n)

    ## compute glm prediction region
    #interval.plugin <- cbind(pred - qnorm(1 - alpha/2) * sepred, 
    #  pred + qnorm(1 - alpha/2) * sepred)
    #sepred <- sd.res * sqrt(1 + 1/n) ## more needed
  #  interval.plugin <- cbind(pred - qnorm(1 - alpha/2) * sepred, 
  #    pred + qnorm(1 - alpha/2) * sepred)
  #}


  ## Get MLEs and plugin interval for gamma gamma distribution 
  shapeMLE <- rateMLE <- 0
  if(family == "Gamma"){
    m1 <- glm(formula, data = data, family = "Gamma")
    #InvFish <- vcov(m1)
    betaMLE <- m1$coefficients
    shapeMLE <- as.numeric(gamma.shape(m1)[1])
    rateMLE <- cbind(1, X) %*% betaMLE * shapeMLE

    #p1 <- predict(m1, type = "response", 
    #  newdata = data.frame(newdata), se.fit = TRUE)
    #pred <- p1$fit
    #sepred <- p1$se.fit
    ## compute glm prediction region 
    #interval.glm <- cbind(pred - qnorm(1 - alpha/2) * sepred, 
    #  pred + qnorm(1 - alpha/2) * sepred)

    ## compute interval prediction region 
    #sepred <- sqrt(1 / shapeMLE * pred^2)
    #interval.plugin <- cbind(pred - qnorm(1 - alpha / 2) * sepred, 
    #  pred + qnorm(1 - alpha / 2) * sepred)
  }


  ## Get MLEs and plugin interval for Inverse Gaussian Distribution 
  if(family == "inverse.gaussian"){
    m1 <- glm(formula, data = data, family = family)
    betaMLE <- m1$coefficients
    #pred <- predict(m1, type = "response", newdata = data.frame(newdata))
    #gx <- 1 / sqrt( cbind(1, X) %*% betaMLE  )
    #scaleMLE <- 1 / mean((Y - gx)^2 / (Y * gx^2) )
    #sepred <- sqrt(pred^3 * scaleMLE)
    
    ## compute glm prediction region 
    #interval.plugin <- cbind(pred - qnorm(1 - alpha/2) * sepred, 
    #  pred + qnorm(1 - alpha/2) * sepred)
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
  index.pred <- findindex(matrix(newdata, ncol = p))
  indices.pred <- sort(unique(index.pred))

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
  newdata.variables <- as.matrix(
    model.matrix(~ ., data.frame(newdata))[, -1])
 

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
      out <- mclapply(1:n.pred, mc.cores = cores, FUN = function(j){
        #x <- newdata.list[[j]]
        x <- newdata.variables[j, ]
        index.j <- which(index.pred[j] == indices.pred)
        #datak <- subkey[[index.j]], ncol = k+1)
        datak <- matrix(subkey[[index.j]], ncol = p + 2)
        Xk <- matrix(datak[, 1:p], ncol = p)
        Yk <- datak[, p+1]
        nk <- nrow(datak)
        rateMLE <- cbind(1, Xk) %*% betaMLE * shapeMLE
        rate.x <- as.numeric(c(1, x) %*% betaMLE * shapeMLE)
        datak.y <- matrix(0, nrow = nrow(datak) + 1, ncol = ncol(datak) - 1)
        colnames(datak.y) <- colnames(data)
        data.y <- matrix(0, nrow = nrow(data) + 1, ncol = ncol(data))
        colnames(data.y) <- colnames(data)
        phatxy <- function(z){
          out <- NULL
          #datak.y[, 1] <- c(Yk, z)
          #datak.y[, -1] <- rbind(Xk, x)
          #datak.y <- as.data.frame(datak.y)
          data.y[, 1] <- c(Y, z)
          data.y[, -1] <- rbind(X, x)
          data.y <- as.data.frame(data.y)          
          m1.y <- glm(formula, data = data.y, family = family)
          betaMLE.y <- m1.y$coefficients

          if(family == "Gamma"){
            ## phat|Ak
            shapeMLE.y <- as.numeric(gamma.shape(m1.y)[1])
            rateMLE.y <- cbind(1, data.y[, -1]) %*% betaMLE.y * shapeMLE.y
            out <- dgamma(data.y[, 1], rate = rateMLE.y, shape = shapeMLE.y)
          }

          if(family == "gaussian"){
            ## phat|Ak
            sd.y <- summary(m1.y)$sigma
            out <- dnorm(data.y[, 1], mean = cbind(1, data.y[, -1]) %*% betaMLE.y, sd = sd.res)
          }

          if(family == "inverse.gaussian"){
            out <- dinvgauss(data.y[, 1], mean = 1 / sqrt(cbind(1, data.y[, -1]) %*% betaMLE.y))
          }

          out
        }

        ## set up range of candidate points used to 
        ## construct the parametric conformal prediction 
        ## region
        lwr <- min(Yk)
        if(min(Yk) < 0) lwr <- 3 * lwr
        if(min(Yk) >= 0) lwr <- (1/3) * lwr
        upr <- max(Yk)
        if(max(Yk) < 0) upr <- (1/3) * upr
        if(max(Yk) >= 0) upr <- 3 * upr

        ## an intial crude approximation of the 
        ## parametric conformal prediction region 
        nk.tilde <- floor(alpha * (nk + 1))
        y.lwr <- lwr
        y.upr <- upr

        if(nk.tilde > 0){

          ## perform a crude crude search 
          crudewidth <- sqrt((upr - lwr)*0.001)
          crude.seq.y <- as.list(seq(from = lwr, to = upr, 
            by = crudewidth))
          crude.search <- lapply(crude.seq.y, 
            FUN = function(y){ 
              foo <- cbind(phatxy(y), c(index.pred, index.j))
              bar <- foo[which(foo[, 2] == index.j), -2]
              baz <- rank(bar)[length(bar)]
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

            if(family == "Gamma") endpt1.lwr <- 0.00001
            precise.seq.y1 <- as.list(seq(from = endpt1.lwr, 
              to = endpt2.lwr, length = floor(crudewidth/0.001)))
            precise.search1 <- lapply(precise.seq.y1, 
              FUN = function(y){
                foo <- cbind(phatxy(y), c(index.pred, index.j))
                bar <- foo[which(foo[, 2] == index.j), -2]
                baz <- rank(bar)[length(bar)]
                baz
              })
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
              FUN = function(y){ 
                foo <- cbind(phatxy(y), c(index.pred, index.j))
                bar <- foo[which(foo[, 2] == index.j), -2]
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

    nonparaconformal <- COPS(newdata)
  }

  out = list(paraconformal = paraconformal, 
    nonparaconformal = nonparaconformal)
  return(out)
}


