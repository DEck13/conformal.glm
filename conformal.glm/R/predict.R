

predict.paraconformal <- function(object, ..., newdata = NULL, alpha = 0.10, 
	cores = 1, bins = NULL, parametric = TRUE, intercept = TRUE, 
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
  }
  newdata <- as.matrix(newdata)
  mf <- model.frame(object)
  Y <- model.response(mf)
  #offs <- model.offset(mf)

  int <- regions(formula = formula, data = data, newdata = newdata, 
  	family = "family", link, alpha = alpha, cores = cores, bins = bins, 
  	intercept = intercept, parametric = parametric, 
    nonparametric = nonparametric)

  paraconformal <- int$paraconformal
  nonparaconformal <- int$nonparaconformal

  out = list(paraconformal = paraconformal, 
    nonparaconformal = nonparaconformal)
  return(out)
}

