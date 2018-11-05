

library(conformalInference)
alpha <- 0.10

## get newdata in form that is accepted by LS conformal 
## and our paramteric conformal implementation
newdata.variables <- as.matrix(
  model.matrix(~ ., data.frame(newdata))[, -1])

## Get regression based conformal prediction region 
## using conformal.pred in the conformalInference 
## package 
funs <- lm.funs(intercept = TRUE)
train.fun <- funs$train.fun
predict.fun <- funs$predict.fun
p1.tibs <- conformal.pred(x = X, y = Y, x0 = newdata.variables, 
  train.fun = train.fun, predict.fun = predict.fun, 
  alpha = alpha, grid.method ="linear",
  num.grid.pts = 999)
LSconformal <- cbind(p1.tibs$lo, p1.tibs$up)

