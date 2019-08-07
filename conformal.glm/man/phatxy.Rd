\name{phatxy}
\title{Conformal Scores for Generalized Linear Regression Models}
\alias{phatxy}
\concept{regression}
\concept{linear regression}
\concept{nonparametric regression}
\concept{Gamma}
\concept{exponential}
\concept{inverse Gaussian}
\concept{exponential family}
\concept{conformal prediction}
\concept{conformal score}
\concept{finite sample validity}
\description{
  Internal function for \code{conformal.glm} package.
  This function computes parametric conformal scores for an augmented 
  dataset which includes the original dataset and a new data point 
  (xnew, ynew).  Functionality is provided for the normal, Gamma, and 
  inverse Gaussian families in the \code{glm} package.  
}
\usage{
  phatxy(ynew, xnew, Yk, Xk, xnew.modmat, 
    data, formula, family, link)
}
\arguments{
  \item{ynew}{candidate response value.}
  \item{xnew}{candidate predictor values.}
  \item{Yk}{the vector of responses which corresponding to predictors 
    that belong to element k of the partition of the predictor space.}
  \item{Xk}{matrix of predictors for which each row is a point within  
    element k of the partition of the predictor space.}
  \item{xnew.modmat}{an internal quantity computed within 
    \code{\link{conformal.glm}}.  A model matrix which includes the 
    predictor value xnew.}
  \item{data}{a data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{regions} is called.}
  \item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.
    See \code{\link[stats]{glm}} and \code{\link[stats]{formula}}
    for description of the R formula mini-language.}
  \item{family}{a character string specifying the family, must be one of
    \code{"gaussian"} (default), \code{"Gamma"}, or 
    \code{"inverse.gaussian"}. May be abbreviated.}
  \item{link}{the function which takes the conditional expectation of 
    the response variable given predictors as its argument and has the 
    linear regression equation as its output.  If missing then the 
    default link function in \code{glm} will be specified.}  
}
\details{
  Internal function for \code{conformal.glm} package.
  This function computes parametric conformal scores for an augmented 
  dataset which includes the original dataset and a new data point 
  (xnew, ynew).  Functionality is provided for the normal, Gamma, and 
  inverse Gaussian families in the \code{glm} package.  
}
\value{
  \item{out}{vector of conformal scores.  The last element of this 
    vector is the conformal score corresponding to the (xnew, ynew) 
    data point.}
}

\seealso{
  \code{\link{conformal.glm}}, \code{glm}
}
\references{
  Eck, D.~J., Crawford, F.~W., and Aronow, P.~M. (2018+)
  Conformal prediction for exponential families and generalized linear models.
  Preprint available on request (email \url{daniel.eck@yale.edu}).

  Lei, J., G'Sell, M., Rinaldo, A., Tibshirani, R., and Wasserman, L. (2016)
  Distribution-Free Predictive Inference for Regression. 
  \url{https://arxiv.org/abs/1604.04173}

  Lei, J. and Wasserman, L. (2014)
  Distribution-Free Prediction Bands for Non-parametric Regression. 
  Journal of the Royal Statistical Society: Series B, 76(1), 71-96.

  Lei, J., Robins, J., and Wasserman, L. (2013)
  Distribution Free Prediction Sets. 
  Journal of the American Statistical Association, 108(501), 278-287.
}
\examples{

}
\keyword{generalized linear regression models}
\keyword{conformal prediction}