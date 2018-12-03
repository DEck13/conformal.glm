\name{local.coverage}
\title{Local coverage calculator}
\alias{local.coverage}
\description{
  Assigns bin labels to data points.
}
\usage{
  local.coverage(region, data, newdata, k, bins = NULL,
    at.data = "TRUE")
}
\arguments{
  \item{region}{a prediction region of interest.}
  \item{data}{a data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{regions} is called. The 
    \code{local.coverage} functionality assumes that the response variable 
    is in the first column.}
  \item{newdata}{an optional matrix, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    new observations for which a prediction is desired. If missing, then 
    prediction regions will be provided at the observed data.}
  \item{k}{desired dimension of the predictor space corresponding to the 
    model of interest.} 
  \item{bins}{an optional argument for specifying the desired number of 
    bins to use along one dimension of the predictor space. If missing, the 
    theoretical large sample optimal bin width is used 
    (width = \code{O(log(n)/n)^(1/(d+1))} where \code{n} is the sample size 
    and \code{d} is the dimension of the main effects).}  
  \item{at.data}{a boolean variable specifying whether or not predictions 
    are made at the original data.}
}
\details{
  This function is a beta version. Functionality is implemented for univariate 
  problems.

  The \code{local.coverage} function estimates the coverage probability 
  of a prediction region locally. We partition the predictor space into 
  equilateral cubes (number of bins is equal to \code{bins^d}) and then 
  estimate the coverage probability of the provided prediction region 
  \code{region} locally for each bin forming the partition.

  When \code{at.data = "FALSE"}, we fit a thrid-order polynomial model 
  to form the lower and upper ends of the prediction region.
}
\value{
  \item{output}{estimated coverage probabilities of the region given by 
    \code{region} for each bin.  If each coverage probability is greater 
    or equal to \code{1 - alpha}, then we say that the region satisties 
    finite sample local validity.}
}

\seealso{
  \code{\link{regions},\link{conformal.glm}}
}