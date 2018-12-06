\name{simulator}
\title{Simulator for comparing the performance of prediction regions}
\alias{simulator}
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
  functionality to construct the parametric conformal prediction region,  
  the nonparametric conformal prediction region, the prediction region 
  formed by the conformalization of residuals, and the highest density 
  region under a specified model. This function generates a dataset 
  and fits a main-effects only model internally.  The output from this 
  function is the estimated local coverage probabilities for each 
  prediction region and the estimated area of each prediction region.
}
\usage{
  simulator(n, alpha = 0.10, beta, bins = NULL, family = "Gamma", 
    link = "inverse", shape = NULL, sd = NULL, confamily = "Gamma",
    parametric = TRUE, nonparametric = FALSE, 
    LS = FALSE, HD = FALSE, cores = 6)
}
\arguments{
  \item{n}{the sample size of an internal dataset.}
  \item{alpha}{the error tolerance desired for the prediction region.  
    The default is set at 0.10.}
  \item{beta}{a regression coefficient vector that is desired for 
    constructing the internal dataset.}
  \item{bins}{an optional argument for specifying the desired number of 
    bins to use along one dimension of the predictor space. If missing, the 
    theoretical large sample optimal bin width is used 
    (width = \code{O(log(n)/n)^(1/(d+1))} where \code{n} is the sample size 
    and \code{d} is the dimension of the main effects).}
  \item{family}{a character string specifying the family used to construct 
    the internal dataset, must be one of
    \code{"Gamma"} (default), \code{"gaussian"}, or 
    \code{"inverse.gaussian"}. May be abbreviated.}
  \item{link}{the function which takes the conditional expectation of 
    the response variable given predictors as its argument and has the 
    linear regression equation as its output.  Only relevant for the 
    \code{"Gamma"} family.  Choices are \code{"inverse"} (default, and 
    consistent with \code{glm} defualt), \code{"identity"}, 
    and \code{"log"}.  If missing then the 
    default link function in \code{glm} will be specified.}  
  \item{shape}{the shape parameter for a \code{"Gamma"} regression.}
  \item{sd}{the standard deviation of the errors used to generate the 
    \code{family = "gaussian"} internal dataset.}
  \item{confamily}{a character string specifying the family used to 
    calculate the parametric conformal prediction region and the 
    highest density region, must be one of
    \code{"Gamma"} (default), \code{"gaussian"}, or 
    \code{"inverse.gaussian"}. May be abbreviated.}
  \item{cores}{the number of cores used to compute the conformal 
    prediction regions. The default is set at 1 core so that the 
    \code{parallel} library is not used. Specifying \code{cores} greater 
    than 1 will call \code{mclapply} package in the \code{parallel} package 
    to speed up the optimization process.}  
  \item{parametric}{a Boolean variable corresponding to whether or not the 
    parametric conformal region is to be computed. The default is set at TRUE.}
  \item{nonparametric}{a Boolean variable corresponding to whether or not the 
    nonparametric conformal region is to be computed. The default is set at 
    TRUE.}
  \item{LS}{a Boolean variable corresponding to whether or not the 
    prediction region formed by the conformalization of residuals is to be 
    computed. The default is set at TRUE.}
  \item{HD}{a Boolean variable corresponding to whether or not the 
    highest density region is to be computed. The default is set at TRUE.}      
}
\details{
  This function first generates \code{n} observations of \code{p} predictors 
  independently from a \code{U(0,1)} where \code{p} is the dimension of 
  \code{beta} minus 1. Then a response variable is generated from the 
  exponential family specified by \code{family}.  For example, we call 
  \code{rgamma(n = n, shape = shape, rate = (cbind(1, x) %*% beta) * shape)} 
  to generate a response variable when \code{family = "Gamma"} and 
  \code{link = "inverse"}. We now have our dataset.

  This function then fits a \code{glm} model to this dataset. 

  If \code{parametric = TRUE}, the parametric conformal prediction region is 
  obtained by calling \code{conformal.glm} on the fitted \code{glm} model.
  If \code{nonparametric = TRUE}, the nonparametric conformal prediction 
  region is obtained by calling \code{conformal.glm} on the fitted \code{glm} 
  model. If \code{LS = TRUE}, the prediction region formed by the 
  conformalization of residuals is obtained by calling the 
  \code{conformal.pred} function from the \code{conformalInference} package 
  (this package needs to be installed for this to work). If \code{HD = TRUE}, 
  then the highest density region is obtained by calling the \code{hdi} 
  function in the \code{HDInterval} package (this package needs to be installed 
  for this to work). All prediction regions are at level 1 - \code{alpha}.  No 
  \code{newdata} is specified, the predicted intervals are for each \code{x} 
  generated within this function.  The parametric and nonparametric conformal 
  prediction regions are computed locally, with \code{bins^(p-1)} specifying 
  the number of equilateral cubes forming the partition of the predictor 
  space. 

  Performance of all prediction regions is assessed by estimated local 
  coverage probabilities corresponding to data within each \code{bin}. 
  The area of each prediction region is also estimated. 
}
\value{
  \item{output.parametric}{estimated local coverage probabilities and area 
    (area is the last element of this vector) of the parametric conformal 
    prediction region.}
  \item{output.nonparametric}{estimated local coverage probabilities and area 
    (area is the last element of this vector) of the nonparametric conformal 
    prediction region.}
  \item{output.LS}{estimated local coverage probabilities and area 
    (area is the last element of this vector) of the prediction region formed 
    by the conformalization of residuals.}
  \item{output.HD}{estimated local coverage probabilities and area 
    (area is the last element of this vector) of the highest density 
    prediction region.}
}
\seealso{
  \code{\link{conformal.glm}}, \code{\link{regions}}, \code{glm}
}
\references{
  Eck, D.~J., Crawford, F.~W., and Aronow, P.~M. (2018+)
  Conformal prediction for exponential families and generalized linear models.
  Preprint available on request (email \url{daniel.eck@yale.edu}).

  Meredith, M. (2018) R package HDInterval.
  \url{https://cran.r-project.org/web/packages/HDInterval/}

  Lei, J., G'Sell, M., Rinaldo, A., Tibshirani, R., and Wasserman, L. (2016)
  Distribution-Free Predictive Inference for Regression. 
  \url{https://arxiv.org/abs/1604.04173}

  Tibshirani, R. (2016) R package conformalInference.
  \url{https://github.com/ryantibs/conformal}

  Lei, J. and Wasserman, L. (2014)
  Distribution-Free Prediction Bands for Non-parametric Regression. 
  Journal of the Royal Statistical Society: Series B, 76(1), 71-96.

  Lei, J., Robins, J., and Wasserman, L. (2013)
  Distribution Free Prediction Sets. 
  Journal of the American Statistical Association, 108(501), 278-287.
}
\examples{
# example of section 2.4 in Geyer (2009)
# data(sports)
# out <- glmdr(cbind(wins, losses) ~ 0 + ., family = "binomial", data = sports)
#summary(out)
}
\keyword{generalized linear regression models}
\keyword{conformal prediction}