\name{find.index}
\title{Support function for the regions function}
\alias{find.index}
\description{
  Assigns bin labels to data points.
}
\usage{
  find.index(mat, wn, k)
}
\arguments{
  \item{mat}{a model matrix consisting of predictors, main effects only.}
  \item{wn}{the bin width corresponding to a partition of the 
    predictor space.  The partitions considered are formed by 
    equilateral cubes.  Therefore \code{wn} is a scalar.}
  \item{k}{desired dimension of the predictor space corresponding to the 
    model of interest.} 
}
\details{
  This function assigns bin labels to data points.  It is used internally in 
  \code{regions.R} and is called on the model matrix of the main effects and 
  the \code{newdata} matrix.
}
\value{
  \item{out}{bin labels}
}

\seealso{
  \code{\link{regions}}
}