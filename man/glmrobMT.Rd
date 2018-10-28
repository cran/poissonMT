\name{glmrobMT}
\alias{glmrobMT}
\title{
  Robust Fitting of Poisson Generalized Linear Models
  using MT robust method
}
\description{
\code{glmrobMT} is used to fit generalized linear models by robust MT
  method. The model is specified by the \code{x} and \code{y}
  components and a description of the error distribution. Currently,
  only implemented for \code{family=poisson}.
}
\usage{
glmrobMT(x, y, weights=NULL, start=NULL, offset=NULL,
  family=poisson(), weights.on.x="none",
  control=glmrobMT.control(), intercept=TRUE, trace.lev=1,
  include.cubinf=TRUE, m.approx=NULL, mprime.approx=NULL, ...)
}
\arguments{
\item{x}{
  design matrix of dimension \eqn{n \times p}{n * p}.
}
\item{y}{
  vector of observations of length \code{n}.
}
\item{weights}{
  an optional vector of weights to be used in the fitting
  process (in addition to the robustness weights computed in
  the fitting process).  
}
\item{start}{
  starting values for the parameters in the linear predictor.
  Note that specifying \code{start} skips the computation of
  the initial estimates, but needs to be \emph{robust} itself.  
}
\item{offset}{
  this can be used to specify an \emph{a priori} known component
  to be included in the linear predictor during fitting.
  At the moment it is not used.
}
\item{family}{
  a description of the error distribution and link function to
  be used in the model. This can be a character string naming a
  family function, a family \code{\link{function}} or the result of a
  call to a family function. See \code{\link{family}} for details of
  family functions. At the moment only \code{poisson} is available.
}
\item{weights.on.x}{
  a character string (can be abbreviated), a \code{\link{function}} or
  \code{\link{list}} (see below), or a numeric vector of length
  \code{n}, specifying how points (potential outliers) in x-space are
  downweighted. If \code{"hat"}, weights on the design of the form
  \eqn{\sqrt{1-h_{ii}}} are used, where \eqn{h_{ii}} are the diagonal
  elements of the hat matrix. If \code{"robCov"}, weights based on
  the robust Mahalanobis distance of the design matrix (intercept
  excluded) are used where the covariance matrix and the centre is
  estimated by \code{\link[MASS]{cov.rob}}
  from the package \pkg{MASS}.\cr
  Similarly,  if \code{"covMcd"}, robust weights are computed using
  \code{\link{covMcd}}.  The default is \code{"none"}.

  If \code{weights.on.x} is a \code{\link{function}}, it is called
  with arguments \code{(X, intercept)} and must return an n-vector of
  non-negative weights.

  If it is a \code{\link{list}}, it must be of length one, and as
  element contain a function much like \code{\link{covMcd}()} or
  \code{\link[MASS]{cov.rob}()} (package \pkg{MASS}), which computes
  multivariate location and \dQuote{scatter} of a data matrix \code{X}.
}
\item{control}{
  a list of parameters for controlling the fitting process.
  See the documentation for \code{\link{glmrobMT.control}} for
  details.
}
\item{intercept}{
  logical indicating if an intercept at the first column of \code{x}
  is present. This information is only used when \code{weights.on.x}
  is not set to \code{none}.
}
\item{trace.lev}{
  logical (or integer) indicating if intermediate results
  should be printed; defaults to \code{0} (the same as \code{FALSE}).
}
\item{include.cubinf}{
  logical, if \code{TRUE} the \code{\link[robcbi]{cubinf}}
  is also used as possibile starting value.
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
\item{\dots}{
  At the moment it is not used.  
}
}
\value{
A list with the following components:
\item{coefficients}{a named vector of coefficients.}
\item{initial}{Initial vector of coefficients.}
\item{family}{the \code{\link{family}} object used.}
\item{residuals}{weighted Pearson residuals.}
\item{fitted.values}{the fitted mean values, obtained by transforming
    the linear predictors by the inverse of the link function.}
\item{linear.predictors}{the linear fit on link scale.}
\item{cov}{the estimated asymptotic covariance matrix of the estimated
    coefficients.}
\item{converged}{logical. Was the IWLS algorithm judged to have converged?}
\item{iter}{the number of iterations used by the influence algorithm.}
\item{cw}{the tuning constant c in Tukey's bisquare psi-function.}
\item{weights.on.x}{how the weights on the design matrix \code{x} were evaluated.}
\item{w.x}{weights used to down-weight observations based on the
    position of the observation in the design space.}
\item{w.r}{robustness weights for each observations.}
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
  \code{\link{poissonMT}}, \code{\link{glmrob}} and \code{\link{cubinf}}
}
\examples{
  data(epilepsy)  
  Efit1 <- glm(Ysum ~ Age10 + Base4*Trt, family=poisson, data=epilepsy)

  x <- model.matrix( ~ Age10 + Base4*Trt, data=epilepsy)
  poissonMTsetwd(tempdir())   
  Efit2 <- glmrobMT(x=x, y=epilepsy$Ysum) 
}
  
\keyword{robust}
\keyword{regression}
\keyword{nonlinear}
