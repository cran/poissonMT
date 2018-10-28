#############################################################
#	poissonMT function
#	Author: C. Agostinelli, M. Valdora and V.J. Yohai
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: 18 August 2018
#	Version: 0.2-2
#	Copyright (C) 2018 C. Agostinelli M. Valdora and V.J. Yohai
#############################################################

poissonMT <- function(x, y, start, weights=NULL, tol=1e-8, maxit=100, m.approx=NULL, mprime.approx=NULL, psi="bisquare", cc=2.3, na.to.zero=TRUE) {
  if (psi!="bisquare")
    stop("Only 'psi==bisquare' since mk.m_rho internal function")
  z <- sqrt(y)
  if (is.null(weights))
    weights <- rep(1,NROW(x))
  good <- weights > 0.0001
  X <- x[good,,drop=FALSE]
  Z <- z[good]
  weights <- weights[good]
  n <- NROW(X)
  p <- NCOL(X)
  difference <- tol+1
  iter <- 0
  if (is.null(m.approx))
    m.approx <- mk.m_rho(cw=cc)
  if (is.null(mprime.approx))
    mprime.approx <- mk.mprime_rho(cw=cc)
  while (difference > tol & iter <= maxit) {
    iter <- iter + 1
    eta <- drop(X%*%start)
    mu <- exp(eta)
    mu[mu > .Machine$double.xmax] <- .Machine$double.xmax    
    w <- mmd(mu, mprime.approx)*mu
    resid <- Z - mm(mu, m.approx)
    wRob <- Mwgt(resid, cc, psi)
    wz <- resid + w*eta
    newstart <- lm.wfit(x=X, y=wz/w, w=w^2*wRob*weights)$coefficients
    newstart[is.na(newstart)] <- 0 ## So that we have a result also when n < p
    difference <- max(abs(start-newstart))
    start <- newstart
  }
  ans <- list()
  start <- lm.wfit(x=X, y=wz/w, w=w^2*wRob*weights)$coefficients
  if (na.to.zero)
    start[is.na(start)] <- 0
  ans$coefficients <- drop(start)
  eta <- drop(x%*%start)
  mu <- exp(eta)
  w <- mmd(mu, mprime.approx)*mu
  ans$fitted.values <- mu
  ans$linear.predictors <- eta 
  ans$residuals <- z - mm(mu, m.approx)
  ans$weights <- w
  ans$w.r <- wRob
  ans$prior.weights <- weights
  ans$converged <- iter < maxit
  ans$iter <- iter
  ans$obj <- sum(Mchi(x=ans$residuals, cc=cc, psi=psi))
  return(ans)
}
