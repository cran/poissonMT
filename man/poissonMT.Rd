\name{poissonMT}
\alias{poissonMT}
\title{
  Robust Fitting of Poisson Generalized Linear Models
  using MT robust method
}
\description{
\code{poissonMT} is used to fit generalized linear models by robust MT
  method. The model is specified by the \code{x} and \code{y}
  components.
}
\usage{
  poissonMT(x, y, start, weights = NULL, tol = 1e-08, maxit = 100,
  m.approx = NULL, mprime.approx = NULL, psi = "bisquare",
  cc = 2.3, na.to.zero = TRUE)
}
\arguments{
\item{x}{
  design matrix of dimension \eqn{n \times p}{n * p}.
}
\item{y}{
  vector of observations of length \code{n}.
}
\item{start}{
 starting values for the parameters in the linear predictor.
}
\item{weights}{
  an optional vector of weights to be used in the fitting
  process (in addition to the robustness weights computed in
  the fitting process).
}
\item{tol}{
  convergence tolerance for the parameter vector.
}
\item{maxit}{
  integer specifying the maximum number of IRWLS iterations.
}
\item{m.approx}{
  a function that return the value, for each linear predictor, that
  makes the estimating equation Fisher consistent. If \code{NULL}
  the default internal function is used.  
}
\item{mprime.approx}{
  a function that return the value, for each linear predictor,
  corresponding to the first derivative of \code{m.approx}.
  If \code{NULL} the default internal function is used. 
}  
\item{psi}{
  the name of the \code{psi} function. At the moment only the
  \code{bisquare} is available.
}
\item{cc}{
  tuning constant c for Tukey's bisquare psi-function.
}
\item{na.to.zero}{
  logical, should the eventual \code{NA} in the coefficients
  be replaced by \code{0}? 
}
}
\value{
A list with the following components
\item{coefficients}{a vector of coefficients.}
\item{fitted.values}{the fitted mean values, obtained by transforming
    the linear predictors by the inverse of the link function.}
\item{linear.predictors}{the linear fit on link scale.}  
\item{residuals}{residuals on the transformed scale.}
\item{weights}{the \emph{working} weights, that is the weights in the
  final iteration of the IWLS fit.}
\item{w.r}{robustness weights for each observations.}    
\item{prior.weights}{the weights initially supplied, a vector of
  \code{1}s if none were.}
\item{converged}{logical. Was the IWLS algorithm judged to have converged?}
\item{iter}{the number of iterations used by the influence algorithm.}
\item{obj}{value of the MT objective function at \code{coefficients}.}
}
\references{
  C. Agostinelli, M. Valdora and V.J Yohai (2018)
  Initial Robust Estimation in Generalized Linear Models
  with a Large Number of Covariates. Submitted.

  M. Valdora and V.J. Yohai (2014)
  Robust estimators for Generalized Linear Models.
  Journal of Statistical Planning and Inference, 146, 31-48.
}
\author{
  Claudio Agostinelli, Marina Valdora and Victor J. Yohai  
}
\seealso{
  \code{\link{glmrobMT}}, \code{\link{glmrob}} and \code{\link{cubinf}}
}
\examples{
  data(epilepsy)
  x <- model.matrix( ~ Age10 + Base4*Trt, data=epilepsy)
  poissonMTsetwd(tempdir())  
  start <- poissonMTinitial(x=x, y=epilepsy$Ysum)$coefficients
  Efit3 <- poissonMT(x=x, y=epilepsy$Ysum, start=start)
}

\keyword{robust}
\keyword{regression}
\keyword{nonlinear}
