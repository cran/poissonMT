#############################################################
#	poissonL2T function
#	Author: C. Agostinelli, M. Valdora and V.J. Yohai
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: 16 August 2018
#	Version: 0.2-1
#	Copyright (C) 2018 C. Agostinelli M. Valdora and V.J. Yohai
#############################################################
  
poissonL2T <- function(x, y, start=NULL, tol=1e-8, maxit=100, m.approx=NULL, mprime.approx=NULL, na.to.zero=TRUE) {
  z <- sqrt(y)
  n <- NROW(x)
  p <- NCOL(x)
  difference <- tol+1
  iter <- 0
  if (is.null(m.approx))
    m.approx <- mk.m_L2()
  if (is.null(mprime.approx))
    mprime.approx <- mk.mprime_L2()
  if (!is.null(start)) {
    eta <- drop(x %*% start)
    mu <- exp(eta)
  } else {
    start <- rep(0, p)
    mu <- y+0.1
    eta <- log(mu)
  }  
  while (difference > tol & iter <= maxit) {
    iter <- iter + 1
    mu[mu > .Machine$double.xmax] <- .Machine$double.xmax
    w <- wL2T(mu, mprime.approx)
    wz <- z - mmL2T(mu, m.approx) + w*eta
    keep <- w > 1e-8
    if (any(keep))
      newstart <- lm.wfit(x=x[keep,], y=(wz/w)[keep], w=(w[keep])^2)$coefficients
##      newstart <- solve(t(x)%*%diag(w^2)%*%x)%*%t(x)%*%diag(w)%*%wz
    else
      newstart <- start <- rep(0, p)
    newstart[is.na(newstart)] <- 0 ## So that we have a result also when n < p
    difference <- max(abs(start-newstart))
####    start <- (start+newstart)/2
    start <- lm.wfit(x=x[keep,], y=(wz/w)[keep], w=(w[keep])^2)$coefficients
    if (na.to.zero)
      start[is.na(start)] <- 0
    eta <- drop(x%*%start)
    mu <- exp(eta)    
  }
  return(drop(start))
}

#############################################################
#	poissonL2Tequ function
#	Author: C. Agostinelli, M. Valdora and V.J. Yohai
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: 16 August 2018
#	Version: 0.1-1
#	Copyright (C) 2018 C. Agostinelli M. Valdora and V.J. Yohai
#############################################################

poissonL2Tequ <- function(x, y, beta, m.approx, mprime.approx) {
  eta <- drop(x%*%beta)
  mu <- exp(eta)
  t((sqrt(y) - mmL2T(mu, m.approx))*mmprimeL2T(mu, mprime.approx)*mu)%*%x
}
