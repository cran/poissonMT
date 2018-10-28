poissonSSinitial <- function (x, y, nsubm, size=ncol(x), cc=2.3,
  psi="bisquare", na.to.zero=TRUE, trace.lev=0) {
  mm <- mk.m_rho(cw=cc)
  weightMTpoisson <- function(beta, x, Y, cc, psi) {
    eta <- drop(x %*% beta)
    resid <- Y - mm(exp(eta))
    goal <- Mchi(x=resid, cc=cc, psi=psi)
    sum(goal)
  }  
  stopifnot(is.matrix(x), (nsubm <- as.integer(nsubm)) >= 1)
  p <- ncol(x)
  n <- nrow(x)
  s2.best <- Inf
  b.best <- rep(NA_real_, p)
  kk <- 0
  for (l in 1:nsubm) {
    if (trace.lev) {
      if (trace.lev > 1) 
        cat(sprintf("%3d:", l))
      else cat(".", if (l%%50 == 0) 
        paste0(" ", l, "\n"))
    }
    i.sub <- sample(x=n, size=size)
    estim0 <- try(glm.fit(x=x, y=y, family=poisson()))
    if (class(estim0)=="try-error") {
      next
    } else {
      estim0 <- estim0$coefficients
      estim0[is.na(estim0)] <- 0
    }
    eta <- drop(x %*% estim0)
    adev <- abs(y * (log(y + (y == 0)) - eta) - (y - exp(eta)))
    if (trace.lev > 1) 
      cat(sprintf(" D=%11.7g ", sum(adev)))
    half <- ceiling(n/2)
    srt.d <- sort(adev, partial = half)
    podador <- adev <= srt.d[half]
    xPod <- x[podador, ]
    yPod <- y[podador]
    fitE <- try(glm.fit(x=xPod, y=yPod, family = poisson()))
    if (class(fitE)=="try-error") {  
      message("glm(.) {inner subsample} error: ", fitE$message)
      if (trace.lev > 1) 
        cat("\n")
    } else {
      betapod <- drop(fitE$coefficients)
      if (all(is.na(betapod))) 
        next
      if (na.to.zero)
        betapod[is.na(betapod)] <- 0         
      kk <- kk + 1
      s2 <- weightMTpoisson(beta=betapod, x=x, Y=sqrt(y), cc=cc, psi=psi)
      if (trace.lev > 1) 
        cat(sprintf("s2=%14.9g", s2))
      if (s2 < s2.best) {
        if (trace.lev > 1) 
          cat(" New best!\n")
          b.best <- betapod
          s2.best <- s2
      } else if (trace.lev > 1) 
        cat("\n")
    }
  }
  list(coefficients = b.best, obj=s2.best, nOksamples = kk)
}
