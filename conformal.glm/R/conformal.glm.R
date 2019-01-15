

conformal.glm <- function(object, ..., newdata = NULL, alpha = 0.10, 
	cores = 1, bins = NULL, parametric = TRUE, intercept = TRUE, 
	nonparametric = FALSE, h = NULL, precision = 0.001){

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
  #mf <- model.frame(object)
  #Y <- model.response(mf)
  #offs <- model.offset(mf)

  stopifnot("glm" %in% class(object))
  stopifnot(family %in% c("Gamma", "gaussian", "inverse.gaussian"))

  int <- regions(formula = formula, data = data, newdata = newdata, 
  	family = family, link = link, alpha = alpha, cores = cores, 
    bins = bins, parametric = parametric, nonparametric = nonparametric, 
    h = h, precision = precision)

  paraconformal <- int$paraconformal
  nonparaconformal <- int$nonparaconformal
  #interval.plugin <- int$interval.plugin

  out = list(paraconformal = paraconformal, 
    nonparaconformal = nonparaconformal)
  return(out)
}

