


library(MASS) 
library(statmod)

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
  sepred <- p1$se.fit

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

