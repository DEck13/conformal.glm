

conformal.glm <- function(object, ..., newdata = NULL, alpha = 0.10, 
	cores = 1, bins = NULL, parametric = TRUE, intercept = TRUE, 
	nonparametric = FALSE, h = NULL, precision = 0.001){


  out <- NULL 
  if(any(attr(object$terms, "dataClasses")[-1] != "factor")){

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

    stopifnot("glm" %in% class(object))
    stopifnot(family %in% c("Gamma", "gaussian", "inverse.gaussian"))

    int <- regions(formula = formula, data = data, 
      newdata = newdata, family = family, link = link, 
      alpha = alpha, cores = cores, bins = bins, 
      parametric = parametric, nonparametric = nonparametric, 
      h = h, precision = precision)

    paraconformal <- int$paraconformal
    nonparaconformal <- int$nonparaconformal
    out = list(paraconformal = paraconformal, 
      nonparaconformal = nonparaconformal)    
  }


  if(any(attr(object$terms, "dataClasses")[-1] == "factor")){

    #############################################
    ## Does not work yet, functionality 
    ## available elsewhere
    #############################################

    ## some important quantities
    object <- update(object, x = TRUE, y = TRUE)
    X <- object$x
    #Y <- object$y

    call <- object$call
    formula <- call$formula
    fam <- object$family
    family <- fam$family
    link <- fam$link
    data <- object$data

    if(is.null(newdata)){
      #newdata <- X[, colnames(X) %in% colnames(data)]
      newdata <- X[, -1]
    }

    stopifnot("glm" %in% class(object))
    stopifnot(family %in% c("Gamma", "gaussian", "inverse.gaussian"))

    int <- regions(formula = formula, data = data, newdata = newdata, 
  	 family = family, link = link, alpha = alpha, cores = cores, 
      bins = bins, parametric = parametric, nonparametric = nonparametric, 
      h = h, precision = precision)

    paraconformal <- int$paraconformal
    nonparaconformal <- int$nonparaconformal

    out = list(paraconformal = paraconformal, 
      nonparaconformal = nonparaconformal)

  }

  return(out)
}

