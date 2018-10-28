\name{poissonMTinitial}
\alias{poissonMTinitial}
\title{
  Initial Robust Estimates based on MT robust method
  for fitting of Poisson Generalized Linear Models
}
\description{
\code{poissonMTinitial} is used to provides a robust initial estimate
  for fit generalized linear models.
  The model is specified by the \code{x} and \code{y} components.
}
\usage{
poissonMTinitial(x, y, stage2 = TRUE, alpha = c(0.025, 0.025),
  tol = 1e-04, cc = 2.3, psi = "bisquare", maxit = 20,
  zero = sqrt(.Machine$double.eps), replace.small = 0.5, start = NULL,
  na.to.zero = TRUE)
}
\arguments{
\item{x}{
  design matrix of dimension \eqn{n \times p}{n * p}.
}
\item{y}{
  vector of observations of length \code{n}.
}
\item{stage2}{
  logical, the second stage should be performed?
}
\item{alpha}{
  quantile orders used in the second stage.
}
\item{tol}{
  convergence tolerance for the parameter vector.
}
\item{cc}{
  tuning constant c for Tukey's bisquare psi-function.
}
\item{psi}{
  the name of the \code{psi} function. At the moment only the
  \code{bisquare} is available.
}
\item{maxit}{
  integer specifying the maximum number of IRWLS iterations.
}
\item{zero}{
  eigenvalues smaller than \code{zero} will be considered exactly equal to 0.
}
\item{replace.small}{
  all the observations \code{y} smaller than \code{replace.small}
  are replaced by \code{replace.small} value.
}
\item{start}{
  eventual starting values, as a reference, for the parameters
  in the linear predictor.
}
\item{na.to.zero}{
  logical, should the eventual \code{NA} in the coefficients
  be replaced by \code{0}?
}
}
\value{
A list with the following components
  \item{coefficients1}{initial value proposed at the end of the first stage.}
  \item{obj1}{value of the MT objective function at \code{coefficients1}.}
  \item{coefficients2a}{initial value proposed at the end of the first part od the second stage.}
  \item{obj2a}{value of the MT objective function at \code{coefficients2a}.}  
  \item{coefficients2b}{initial value proposed at the end of the second part od the second stage.}
  \item{obj2b}{value of the MT objective function at \code{coefficients2b}.}
  \item{coefficients}{initial value proposed.}
  \item{obj}{value of the MT objective function at \code{coefficients}.}    
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
  \code{\link{poissonMT}} and \code{\link{poissonL2T}}
}
  
\examples{
  data(epilepsy)
  x <- model.matrix( ~ Age10 + Base4*Trt, data=epilepsy)
  poissonMTsetwd(tempdir())
  start <- poissonMTinitial(x=x, y=epilepsy$Ysum)$coefficients
  start
}

\keyword{robust}
\keyword{regression}
\keyword{nonlinear}
