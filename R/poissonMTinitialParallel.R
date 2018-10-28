#############################################################
#	poissonMTinitialParallel function
#	Author: C. Agostinelli, M. Valdora and V.J. Yohai
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: August, 16, 2018
#	Version: 0.2-2
#	Copyright (C) 2018 C. Agostinelli M. Valdora and V.J. Yohai
#############################################################

poissonMTinitialParallel <- function(x, y, stage2=TRUE,
  alpha=c(0.025,0.025), tol=1e-4, cc=2.3, psi="bisquare",
  maxit=20, zero=sqrt(.Machine$double.eps), replace.small=0.5, start=NULL,
  na.to.zero=TRUE, parallel = c("no", "multicore", "snow"), ncpus=1, cl = NULL) {

  if (missing(parallel)) 
    parallel <- "no"
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
    loadNamespace("parallel")
  }
  mkcl <- FALSE
  if (have_snow & is.null(cl)) {
    cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
    mkcl <- TRUE
  }
  n <- NROW(x)
  p <- ncol(x)
  y[y < replace.small] <- replace.small   
  Y <- sqrt(y) ### the data are transformed ones for all!!!
  result <- list()
  mm <- mk.m_rho(cw=cc)
  m.approx <- mk.m_L2()
  mprime.approx <- mk.mprime_L2()  
  
  weightMTpoisson <- function(beta,x,Y,cc,psi) {
    eta <- drop(x %*% beta)
    resid <- Y - mm(exp(eta))
    goal <- Mchi(x=resid, cc=cc, psi=psi)
    sum(goal)
  }

  check <- function(x) {
    !(any(is.na(x) | is.nan(x) | is.infinite(x) | any(x > sqrt(.Machine$double.xmax))))
  }
  ## STEP 1
  start0 <- poissonL2T(x=x, y=y, tol=tol, maxit=maxit, start=start, m.approx=m.approx, mprime.approx=mprime.approx)  
  eta <- drop(x%*%start0)
  mu <- exp(eta)
  mu[mu > .Machine$double.xmax] <- .Machine$double.xmax    
  resid0 <- Y - mmL2T(mu, m.approx)  
  w <- wL2T(mu, mprime.approx)
  wz <- resid0 + w*eta
  keep <- w > sqrt(.Machine$double.eps)
  zz <- (wz/w)[keep]
  xx <- x[keep,]
  yy <- y[keep]
  YY <- Y[keep]
  ww <- (w^2)[keep]
  resid0 <- resid0[keep]  
  nz <- length(zz)
  fndrop1 <- function(i, xx, yy, ww, m.approx, mprime.approx) {    
    beta1 <- lm.wfit(x=xx[-i,], y=zz[-i], w=ww[-i])$coefficients
    beta1[is.na(beta1)] <- 0
    mu1 <- exp(drop(xx%*%beta1))
    mu1[mu1 > .Machine$double.xmax] <- .Machine$double.xmax      
    w1 <- wL2T(mu1, mprime.approx)
    resid1 <- YY - mmL2T(mu1, m.approx)
    return(resid1)
  }
  resid1 <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      as.matrix(as.data.frame(parallel::mclapply(X=1:n, FUN=fndrop1, xx=xx, yy=yy, ww=ww, m.approx=m.approx, mprime.approx=mprime.approx, mc.cores = ncpus)))
    } else if (have_snow) {
      parallel::clusterMap(cl=cl, fun=fndrop1, 1:n, MoreArgs = list(xx=xx, yy=yy, ww=ww, m.approx=m.approx, mprime.approx=mprime.approx), RECYCLE=FALSE, SIMPLIFY=TRUE)
    }
  } else
    as.matrix(as.data.frame(lapply(X=1:n, FUN=fndrop1, xx=xx, yy=yy, ww=ww, m.approx=m.approx, mprime.approx=mprime.approx)))
  R <- t(resid0 - resid1)
  posin <- which(apply(R,1,check))
  k <- length(posin)
  if (k < nz)
    R <- R[posin,,drop=FALSE]
  M <- .Fortran("outerR",
    as.matrix(R),
    M=mat.or.vec(nz,nz),
    as.integer(k),            
    as.integer(nz),
    PACKAGE="poissonMT"
  )$M
  M[M > .Machine$double.xmax/10] <- .Machine$double.xmax/10  
  eM <- eigen(M)
  pos <- eM$values > zero
  vectors <- eM$vectors[,pos]
  Z <- R%*%vectors
  Z <- Z[,apply(Z, 2, sd) > sqrt(.Machine$double.eps), drop=FALSE]    
  zp <- ncol(Z)
  if (zp > 0) {  
    ## Construct A
    xx <- xx[posin,,drop=FALSE]
    yy <- yy[posin]
    ww <- ww[posin]
    fnA <- function(j, Z, xx, yy, ww=ww) {
      mZj <- median(Z[,j])
      A <- matrix(0, nrow=ncol(xx), ncol=3)
      subset1 <- Z[,j] <= mZj
      A[,1] <- lm.wfit(x=xx[subset1,,drop=FALSE], y=zz[subset1], w=ww[subset1])$coefficients
      A[,2] <- lm.wfit(x=xx[!subset1,,drop=FALSE], y=zz[!subset1], w=ww[!subset1])$coefficients 
      aZj <- abs(Z[,j])
      subset1 <- aZj <= median(aZj)
      A[,3] <- lm.wfit(x=xx[subset1,,drop=FALSE], y=zz[subset1], w=ww[subset1])$coefficients
      return(A)
    }
    A <- if (ncpus > 1L && (have_mc || have_snow)) {
      if (have_mc) {
        as.matrix(as.data.frame(parallel::mclapply(X=1:zp, FUN=fnA, Z=Z, xx=xx, yy=yy, ww=ww, mc.cores = ncpus)))
      } else if (have_snow) {
        parallel::clusterMap(cl=cl, fun=fnA, 1:zp, MoreArgs = list(Z=Z, xx=xx, yy=yy, ww=ww), RECYCLE=FALSE, SIMPLIFY=TRUE)
      }
    } else
      as.matrix(as.data.frame(lapply(X=1:zp, FUN=fnA, Z=Z, xx=xx, yy=yy, ww=ww)))
    A <- cbind(A, start0)
    A[is.na(A)] <- 0  
    suma <- rep(0, ncol(A))
    for (j in 1:ncol(A)) {
      suma[j] <- weightMTpoisson(A[,j],x=x,Y=Y,cc=cc,psi=psi)
    }
    coefficients <- A[,which.min(suma)]
    if (stage2) {
      ## Stage 2
      eta <- drop(x%*%coefficients)
      mu <- exp(eta)
      mu[is.infinite(mu)] <- .Machine$double.xmax/10
      w <- wL2T(mu, mprime.approx)
      wz <- Y - mmL2T(mu, m.approx) + w*eta   
      q1 <- qpois(p=alpha[1]/2,mu)
      q2 <- qpois(p=1-alpha[1]/2,mu)
      inside1 <- (y >= q1 & y <= q2)
      keep <- (w > 1e-8) & inside1
      if (sum(keep) > p) {
        z <- (wz/w)[keep]
        xx <- x[keep,]
        w <- (w^2)[keep]
        coefficients2a <- lm.wfit(x=xx, y=z, w=w)$coefficients
        eta <- drop(x%*%coefficients2a)
        mu <- exp(eta)
        mu[mu > .Machine$double.xmax/10] <- .Machine$double.xmax/10    
        w <- wL2T(mu, mprime.approx)
        wz <- Y - mmL2T(mu, m.approx) + w*eta    
        q1 <- qpois(p=alpha[2]/2,mu)
        q2 <- qpois(p=1-alpha[2]/2,mu)
        inside2 <- y >= q1 & y <= q2
        keep <- (w > sqrt(.Machine$double.eps)) & (inside1 | inside2)
        if (sum(keep) > p) {
          z <- (wz/w)[keep]
          xx <- x[keep,]
          w <- (w^2)[keep]
          coefficients2b <- lm.wfit(x=xx, y=z, w=w)$coefficients
        } else {
          coefficients2b <- coefficients2a
          warnings("Stage 2 (b) skipped, since not enough observations in the interval")
        }
      } else {
        coefficients2a <- coefficients2b <- coefficients
        warnings("Stage 2 (a and b) skipped, since not enough observations in the interval")
      }          
    } else {
      coefficients2a <- coefficients2b <- rep(NA,p)
    }
  } else {
    coefficients <- start0
    if (stage2) {
      coefficients2a <- coefficients2b <- coefficients
    } else {
      coefficients2a <- coefficients2b <- rep(NA,p)
    }
  }
  if (have_snow & mkcl)
    parallel::stopCluster(cl)   
  result$coefficients1 <- coefficients
  if (na.to.zero)  
    result$coefficients1[is.na(result$coefficients1)] <- 0
  result$obj1 <- weightMTpoisson(beta=result$coefficients1,x=x,Y=Y,cc=cc,psi=psi)
  result$coefficients2a <- coefficients2a
  if (na.to.zero)
    result$coefficients2a[is.na(result$coefficients2a)] <- 0
  result$obj2a <- weightMTpoisson(beta=result$coefficients2a,x=x,Y=Y,cc=cc,psi=psi)
  result$coefficients2b <- coefficients2b  
  if (na.to.zero)
    result$coefficients2b[is.na(result$coefficients2b)] <- 0
  result$obj2b <- weightMTpoisson(beta=result$coefficients2b,x=x,Y=Y,cc=cc,psi=psi)
  obj <- c(result$obj1, result$obj2a, result$obj2b)
  result$coefficients <- rbind(result$coefficients1, result$coefficients2a, result$coefficients2b)[which.min(obj),]
  result$obj <- obj[which.min(obj)]  
  return(result)
}
