\name{poissonL2T}
\alias{poissonL2T}
\title{
  Fitting of Poisson Generalized Linear Models
  using MT method with L2 rho function
}
\description{
\code{poissonL2T} is used to fit generalized linear models by MT
  method with L2 rho function. The model is specified by the
  \code{x} and \code{y} components. Since the L2 rho function is
  used the method is not robust. 
}
\usage{
  poissonL2T(x, y, start = NULL, tol = 1e-08, maxit = 100,
  m.approx = NULL, mprime.approx = NULL, na.to.zero = TRUE)
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
\item{na.to.zero}{
  logical, should the eventual \code{NA} in the coefficients
  be replaced by \code{0}? 
}  
}
\value{
  A vector with the estimated coefficients.
}
\references{
  C. Agostinelli, M. Valdora and V.J Yohai (2018)
  Initial Robust Estimation in Generalized Linear Models
  with a Large Number of Covariates. Submitted.

  M. Valdora and V.J. Yohai (2014)
  Robust estimators for generalized linear models.
  Journal of Statistical Planning and Inference, 146, 31-48.
}
\author{
  Claudio Agostinelli, Marina Valdora and Victor J. Yohai
}
\seealso{
  \code{\link{poissonMT}}
}
\examples{
  data(epilepsy)
  x <- model.matrix( ~ Age10 + Base4*Trt, data=epilepsy)
  poissonMTsetwd(tempdir())
  Efit4 <- poissonL2T(x=x, y=epilepsy$Ysum)
}
  
\keyword{regression}
\keyword{nonlinear}
