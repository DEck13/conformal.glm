


library(MASS)
library(statmod)
library(parallel)
library(conformal.glm)

## Works but is slow
data(insurance)
m1 <- glm(Y ~ ., data = insurance, family = "Gamma")
gout <- conformal.glm(m1, bins = 3, cores = 6)
gout

## Works but is slow
gout2 <- conformal.glm(m1, bins = 3, parametric = FALSE,  
	nonparametric = TRUE, cores = 6)
gout2

