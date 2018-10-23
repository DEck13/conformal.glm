

conformalprediction <- function(object, ..., newdata = NULL, alpha = 0.10, 
	cores = 6, bins = NULL, parametric = TRUE, LS = FALSE, intercept = TRUE, 
	nonparametric = FALSE){

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
  	family = family, link, alpha = alpha, cores = cores, bins = bins, 
  	intercept = intercept, parametric = parametric, 
    LS = LS, nonparametric = nonparametric)

  paraconformal <- int$paraconformal
  nonparaconformal <- int$nonparaconformal
  LSconformal <- int$LSconformal
  interval.plugin <- int$interval.plugin

  out = list(paraconformal = paraconformal, 
    nonparaconformal = nonparaconformal,
    LSconformal = LSconformal, 
    interval.plugin = interval.plugin)
  return(out)
}

