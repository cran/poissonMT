### Taken from robustbase version 0.93.1 and modified
glmrobMT <- function (x, y, weights=NULL, start=NULL, offset=NULL, 
  family=poisson(), weights.on.x="none", control=glmrobMT.control(), 
  intercept=TRUE, trace.lev=1, include.cubinf=TRUE,
  m.approx=NULL, mprime.approx=NULL, ...) {
  stopifnot(is.numeric(cw <- control$cw), cw > 0)  
  if (family$family != "poisson") 
    stop("Currently, only family 'poisson' is supported for the \"MT\" estimator")
  n <- nrow(x)
  p <- ncol(x)
  if (is.null(weights)) 
    weights <- rep.int(1, n)
  else if (any(weights <= 0)) 
    stop("All weights must be positive")
  if (!is.null(offset)) 
    .NotYetUsed("offset")
  linkinv <- family$linkinv
  variance <- family$variance
  ni <- as.vector(weights)
  sni <- sqrt(ni)
  w <- robXweights(weights.on.x, x, intercept = intercept)
  if (is.null(m.approx))
    m.approx <- mk.m_rho(cw=cw)
  if (is.null(mprime.approx))
    mprime.approx <- mk.mprime_rho(cw=cw)  
  if (is.null(start)) {
    if (trace.lev)
      cat("Computing initial estimate\n")
    start <- poissonMTinitial(x=sqrt(w*weights)*x, y=y*sqrt(w*weights), stage2=TRUE, alpha=c(0.05,0.025,0.025), tol=1e-4, cc=cw, psi="bisquare", maxit=20, zero=sqrt(.Machine$double.eps), replace.small=0.5, start=NULL, na.to.zero=TRUE)$coefficients
    weightMTpoisson <- function(beta,x,y) {
      eta <- drop(x %*% beta)
      resid <- sqrt(y) - mm(exp(eta), m.approx)
      goal <- Mchi(x=resid, cc=cw, psi="bisquare")
      sum(goal)
    }
    if (include.cubinf) {
      ctrl <- cubinf.control()
      start2 <- cubinf(x=sqrt(w)*x,y=y*sqrt(w),weights=weights,family=poisson(), control=ctrl)$coefficients 
      obj1 <- weightMTpoisson(start,x=sqrt(w*weights)*x,y=y*sqrt(w*weights))
      obj2 <- weightMTpoisson(start,x=sqrt(w*weights)*x,y=y*sqrt(w*weights))
      start <- if(obj1 <= obj2) start else start2
    }  
  } else {
    if (!is.numeric(start) || length(start) != p) 
      stop(gettextf("'start' must be an initial estimate of beta, of length %d", p), domain = NA)
  }
  MTres <- poissonMT(x=sqrt(w)*x, y=sqrt(w)*y, start=start, weights=weights, tol=1e-8, maxit=control$maxit, m.approx=m.approx, mprime.approx=mprime.approx, psi="bisquare", cc=cw, na.to.zero=TRUE)
  cov <- covasin(x=x, y=y, beta=MTres$coefficients, cw=cw, m.approx=m.approx, w=w*weights)
  eta <- drop(x %*% MTres$coefficients)
  mu <- linkinv(eta)
  Vmu <- variance(mu)
  if (any(is.na(Vmu))) stop("NAs in V(mu)")
  if (any(Vmu == 0)) stop("0s in V(mu)")
  sVF <- sqrt(Vmu)
  residP <- (y - mu) * sni/sVF
  w.r <- Mwgt(sqrt(y) - mm(exp(eta), m.approx), cw, psi="bisquare")
  names(mu) <- names(eta) <- names(residP)
  names(MTres$coefficients) <- names(start) <- nmB <- colnames(x)
  result <- list(coefficients=MTres$coefficients, initial=start,
    family=poisson(), residuals=residP, fitted.values=mu,
    linear.predictors=eta, cov=cov,
    converged=MTres$converged, iter=MTres$iter, 
    cw=cw, weights.on.x=weights.on.x, w.x=w, w.r=w.r)
  return(result)
}
