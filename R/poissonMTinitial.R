#############################################################
#	poissonMTinitial function
#	Author: C. Agostinelli, M. Valdora and V.J. Yohai
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: 21 August 2018
#	Version: 0.4-3
#	Copyright (C) 2018 C. Agostinelli, M. Valdora and V.J. Yohai
#############################################################

poissonMTinitial <- function(x, y, stage2=TRUE,
  alpha=c(0.025,0.025), tol=1e-4, cc=2.3, psi="bisquare",
  maxit=20, zero=sqrt(.Machine$double.eps), replace.small=0.5, start=NULL,
  na.to.zero=TRUE) {
  n <- NROW(x)
  p <- ncol(x)
  nn <- names(x)
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
  ## STEP 1
  start0 <- poissonL2T(x=x, y=y, tol=tol, maxit=maxit, start=start, m.approx=m.approx, mprime.approx=mprime.approx)
  eta <- drop(x%*%start0)
  mu <- exp(eta)
  mu[mu > .Machine$double.xmax/10] <- .Machine$double.xmax/10
  resid0 <- Y - mmL2T(mu, m.approx)  
  w <- wL2T(mu, mprime.approx)
  wz <- resid0 + w*eta
  keep <- w > sqrt(.Machine$double.eps)
  zz <- (wz/w)[keep]
  xx <- x[keep,,drop=FALSE]
  yy <- y[keep]
  YY <- Y[keep]
  ww <- (w^2)[keep]
  resid0 <- resid0[keep]
  nz <- length(zz)
  posin <- rep(0, nz)
  R <- matrix(NA, nrow=nz, ncol=nz)
  k <- 0
  for (i in 1:nz) {
    beta1 <- lm.wfit(x=xx[-i,], y=zz[-i], w=ww[-i])$coefficients
    beta1[is.na(beta1)] <- 0
    mu1 <- exp(drop(xx%*%beta1))
    mu1[mu1 > .Machine$double.xmax/10] <- .Machine$double.xmax/10 
    w1 <- wL2T(mu1, mprime.approx)
    resid1 <- YY - mmL2T(mu1, m.approx)
    if (!(any(is.na(resid1) | is.nan(resid1) | is.infinite(resid1) | any(resid1 > sqrt(.Machine$double.xmax/10))))) {    
      k <- k + 1
      R[k,] <- resid0 - resid1
      posin[k] <- i
    }
  }
  if (k < nz)
    R <- R[1:k,,drop=FALSE]
  posin <- posin[1:k]
  if (k < p) {
    result$coefficients1 <- rep(NA,p)
    result$opt1 <- NA
    result$coefficients2a <- rep(NA,p)
    result$opt2a <- NA
    result$coefficients2b <- rep(NA,p)
    result$opt2b <- NA    
    return(result)
  }
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
    nA <- 3*zp+1
    A <- matrix(0, nrow=nA, ncol=p)
    for (j in 1:zp) {
      mZj <- median(Z[,j])  
      subset1 <- posin[Z[,j] <= mZj]
      A[3*(j-1)+1,] <- lm.wfit(x=xx[subset1,,drop=FALSE], y=zz[subset1], w=ww[subset1])$coefficients
      subset1 <- posin[Z[,j] >= mZj]
      A[3*(j-1)+2,] <- lm.wfit(x=xx[subset1,,drop=FALSE], y=zz[subset1], w=ww[subset1])$coefficients 
      aZj <- abs(Z[,j])
      subset1 <- posin[aZj <= median(aZj, na.rm=TRUE)]
      A[3*(j-1)+3,] <- lm.wfit(x=xx[subset1,,drop=FALSE], y=zz[subset1], w=ww[subset1])$coefficients     
    }
    A[nA,] <- start0
    A[is.na(A)] <- 0
    suma <- rep(0, nrow(A))
    for (j in 1:nrow(A)) {
      suma[j] <- weightMTpoisson(A[j,],x=x,Y=Y,cc=cc,psi=psi)
    }
    coefficients <- A[which.min(suma),]
    if (stage2) {
      ## Stage 2
      eta <- drop(x%*%coefficients)
      mu <- exp(eta)
      mu[mu > .Machine$double.xmax/10] <- .Machine$double.xmax/10
      w <- wL2T(mu, mprime.approx)
      wz <- Y - mmL2T(mu, m.approx) + w*eta   
      q1 <- qpois(p=alpha[1]/2,mu)
      q2 <- qpois(p=1-alpha[1]/2,mu)
      inside1 <- (y >= q1 & y <= q2)
      keep <- (w > sqrt(.Machine$double.eps)) & inside1
      if (sum(keep) > p) {
        z <- (wz/w)[keep]
        xx <- x[keep,]
        w <- (w^2)[keep]
        coefficients2a <- lm.wfit(x=xx, y=z, w=w)$coefficients
        coefficients2a[is.na(coefficients2a)] <- 0
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
          coefficients2b[is.na(coefficients2b)] <- 0          
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
  names(result$coefficients) <- names(result$coefficients1) <- names(result$coefficients2a) <- names(result$coefficients2b) <- nn
  return(result)
}
