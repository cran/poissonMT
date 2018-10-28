\name{poissonSSinitial}
\alias{poissonSSinitial}
\title{
  Initial Robust Estimates based on SubSampling method
  for fitting of Poisson Generalized Linear Models
}
\description{
\code{poissonSSinitial} is used to provides a robust initial estimate
  for fit generalized linear models.
  The model is specified by the \code{x} and \code{y} components.
}
\usage{
poissonSSinitial(x, y, nsubm, size = ncol(x), cc = 2.3,
  psi = "bisquare", na.to.zero = TRUE, trace.lev = 0)
}
\arguments{
\item{x}{
  design matrix of dimension \eqn{n \times p}{n * p}.
}
\item{y}{
  vector of observations of length \code{n}.
}
\item{nsubm}{
  the number of subsamples to take for finding an initial estimate.
}
\item{size}{
  size of the subsamples.
}  
\item{cc}{
  tuning constant c for Tukey's bisquare psi-function.
}
\item{psi}{
  the name of the \code{psi} function. At the moment only the
  \code{bisquare} is available.
}
\item{na.to.zero}{
  logical, should the eventual \code{NA} in the coefficients
  be replaced by \code{0}?
}
\item{trace.lev}{
  logical (or integer) indicating if intermediate results
  should be printed; defaults to \code{0} (the same as \code{FALSE}).  
}  
}
\value{
A list with the following components
  \item{coefficients}{initial value proposed.}
  \item{obj}{value of the MT objective function at \code{coefficients}.}
  \item{nOksamples}{number of subsamples without errors.}
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
  \code{\link{poissonMTinitial}}
}
  
\examples{
  data(epilepsy)
  x <- model.matrix( ~ Age10 + Base4*Trt, data=epilepsy)
  start <- poissonSSinitial(x=x, y=epilepsy$Ysum, nsubm=100)$coefficients
  start
}

\keyword{robust}
\keyword{regression}
\keyword{nonlinear}
