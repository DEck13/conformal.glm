\name{conformal.glm}
\title{Prediction Regions for Generalized Linear Regression Models}
\alias{conformal.glm}
\concept{regression}
\concept{linear regression}
\concept{nonparametric regression}
\concept{Gamma}
\concept{exponential}
\concept{inverse Gaussian}
\concept{exponential family}
\concept{prediction regions}
\concept{conformal prediction}
\concept{finite sample validity}
\description{
  Compute and compare prediciton regions for the normal, Gamma, and 
  inverse Gaussian families in the \code{glm} package.  There is 
  functionality to construct the parametric conformal prediction 
  region and the nonparametric conformal prediction region. 
}
\usage{
  conformal.glm(object, ..., newdata = NULL, alpha = 0.10, 
  	cores = 1, bins = NULL, parametric = TRUE, intercept = TRUE, 
  	nonparametric = FALSE, h = NULL)
}
\arguments{
  \item{object}{an object of class \code{"glm"}.}
  \item{\dots}{further arguments passed to or from other methods.}
  \item{newdata}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    new observations for which a prediction is desired. If missing, then 
    prediction regions will be provided at the observed data.}  
  \item{alpha}{the error tolerance desired for the prediction region.  
    The default is set at 0.10.}
  \item{cores}{the number of cores used to compute the conformal 
    prediction regions. The default is set at 1 core so that the 
    \code{parallel} library is not used. Specifying \code{cores} greater 
    than 1 will call \code{mclapply} package in the \code{parallel} package 
    to speed up the optimization process.}  
  \item{bins}{an optional argument for specifying the desired number of 
    bins to use along one dimension of the predictor space. If missing, the 
    theoretical large sample optimal bin width is used 
    (width = \code{O(log(n)/n)^(1/(d+1))} where \code{n} is the sample size 
    and \code{d} is the dimension of the main effects).}
  \item{parametric}{a Boolean variable corresponding to whether or not the 
    parametric conformal region is to be computed. The default is set at TRUE.}
  \item{intercept}{a Boolean variable corresponding to whether or not the 
    intercept is included in the regression equation. This is only relevant 
    for the computation of the prediction region by conformalization of 
    residuals (when \code{LS = TRUE}). The default is set at TRUE.}    
  \item{nonparametric}{a Boolean variable corresponding to whether or not the 
    nonparametric conformal region is to be computed. The default is set at 
    TRUE.}
  \item{h}{The bandwidth corresponding to the estimated density to form the 
    nonparametric conformal prediction region.  The default value is 
    \code{1/bins}.}
}
\details{
  This function calls on the \code{\link{regions}} function to compute all of 
  the prediction regions outlined in the description. This function is easier 
  to use than the \code{\link{regions}} function since it can be called 
  directly on an object of class \code{glm}.
}
\value{
  \code{regions} has functionality for the parametric conformal prediction 
  region and the nonparametric conformal prediction region. 

  \item{paraconformal}{The parametric conformal prediction region which is 
    returned when \code{parametric = TRUE}.}
  \item{nonparaconformal}{The nonparametric conformal prediction region which 
    is returned when \code{nonparametric = TRUE}.}
}
\seealso{
  \code{\link{regions}}, \code{glm}
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

