


library(MASS)
library(statmod)
library(conformalInference)
library(parallel)
library(paraconformal)

data(insurance)
m1 <- glm(Y ~ ., data = insurance, family = "Gamma")
gout <- conformalprediction(m1, bins = 3)
gout

