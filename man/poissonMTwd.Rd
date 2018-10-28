\name{poissonMTwd}
\alias{poissonMTsetwd}
\alias{poissonMTgetwd}  
\title{
  Set and get working directory for the package.
}
\description{
  Set and get working directory for the package. This working
  directory is used to store some objects needed during the
  estimation process.
}
\usage{
  poissonMTsetwd(path)

  poissonMTgetwd()
}
\arguments{
\item{path}{
  A valid path name for your OS where you have write/read permission.
}
}

\details{
  \code{poissonMTsetwd} is just a wrapper for \code{options("poissonMT:wd" = path)} while \code{poissonMTsetwd} is just a wrapper for \code{getOption("poissonMT:wd", NULL)}.
}

\value{
  \code{poissonMTsetwd} does not return any object, \code{poissonMTsetwd} return a path or \code{NULL}.  
}
  
\author{
  Claudio Agostinelli, Marina Valdora and Victor J. Yohai  
}
  
\seealso{
  \code{\link{glmrobMT}}, \code{\link{glmrob}} and \code{\link{cubinf}}
}
\examples{
  poissonMTsetwd(tempdir())
  poissonMTgetwd()
}

\keyword{robust}
\keyword{regression}
\keyword{nonlinear}
