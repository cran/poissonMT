#############################################################
#	functions for poissonL2T
#	Author: C. Agostinelli, M. Valdora and V.J. Yohai
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: 16 October 2018
#	Version: 0.1-2
#	Copyright (C) 2018 C. Agostinelli, M. Valdora and V.J. Yohai
#############################################################

mmL2Tint <- function(lambda, tt=sqrt, p=0.9999999, lam.max=100) {
  internal <- function(lambda, tt, p) {
    if (length(lambda)) {
      q <- sapply(qpois(p, lambda), function(x) 0:x, simplify = FALSE)
      mm <- sapply(1:length(lambda), function(i) drop(tt(q[[i]])%*%dpois(q[[i]], lambda[i])))
    } else {
      mm <- vector(length=0)  
    }
    return(mm)
  }
  mm <- rep(0,length(lambda))
  pos <- lambda < lam.max
  mm[pos] <- internal(lambda[pos], tt, p)
  mm[!pos] <- sqrt(lambda[!pos])
  return(mm)
}

mk.m_L2 <- function(lambda = c(seq(0,2.9, by=0.1), seq(3,200)),
  reltol = sqrt(.Machine$double.eps),
  sFile = NULL,
  path = poissonMTgetwd(),
  recompute = getOption("poissonMT:m_rho_recompute", FALSE)) {
  if (is.null(path)) {
    stop("please set up the working directory for this package by calling 'poissonMTsetwd'")
  }
  if (recompute) {
    useFile <- FALSE
    if (is.null(sFile))
      sFile <- file.path(path, "L2Tes.rda")
  } else {
    useFile <- !is.null(sFile)
    if (useFile) {
      sFile <- file.path(path, sFile)
      if (!testPathForOutput(sFile, overwrite=FALSE))
        stop("the 'path' and/or 'sFile' provided are not correct or no permission to write in this file/path")
        load(sFile) #-> 'm.approx'
    } else {
      MML2Tenv <- new.env()
      suppressWarnings(data("L2Tes", envir=MML2Tenv))
      if (exists("m.approx", envir=MML2Tenv, inherits=FALSE)) {
        m.approx <- with(MML2Tenv, m.approx)
        useFile <- TRUE
      } else {
        useFile <- FALSE
        if (is.null(sFile))
          sFile <- file.path(path, "L2Tes.rda")
        else
          sFile <- file.path(path, sFile)  
      }
    }
  }
  if (!useFile) {
    mm.la <- mmL2Tint(lambda, lam.max=max(lambda))
##    m.approx <- approxfun(lambda, mm.la)    
    m.approx <- splinefun(lambda, mm.la, method = "monoH.FC")
    e <- environment(m.approx)
    assign("lambda.max", max(lambda), envir=e)
    if (!testPathForOutput(sFile, overwrite=TRUE))
      stop("the 'path' and/or 'sFile' provided are not correct or no permission to write in this file/path")
    save(m.approx, file = sFile)
  }
  return(m.approx)
}

mmL2T <- function(lam, m.approx) {
  la.max <- environment(m.approx)$lambda.max
  z <- ((m <- lam) <= la.max)
  m[z] <- m.approx(lam[z])
  if(any(i <- !z)) m[i] <- sqrt(lam[i])
  m
}

mmprimeL2Tint <- function(lambda, tt=sqrt, p=0.9999999, lam.max=200) {
  internal <- function(lambda, tt, p) {
    if (length(lambda)) {      
      q <- sapply(qpois(p, lambda), function(x) 0:x, simplify = FALSE)
      mm <- sapply(1:length(lambda), function(i) drop((tt(q[[i]]+1)-tt(q[[i]]))%*%dpois(q[[i]], lambda[i])))
    } else {
      mm <- vector(length=0)  
    }
    return(mm)
  }
  mm <- rep(0,length(lambda))
  pos <- lambda < lam.max
  mm[pos] <- internal(lambda[pos], tt, p)
  mm[!pos] <- 1/(2*sqrt(lambda[!pos]))
###  mm <- ifelse (lambda < 100, internal(lambda, tt, p), 1/(2*sqrt(lambda))) 
  return(mm)  
}

mk.mprime_L2 <- function(lambda = c(seq(0,2.9, by=0.1), seq(3,200)),
  reltol = sqrt(.Machine$double.eps),
  sFile = NULL,
  path = poissonMTgetwd(),
  recompute  = getOption("poissonMT:m_rho_recompute", FALSE)) {
  if (is.null(path)) {
    stop("please set up the working directory for this package by calling 'poissonMTsetwd'")
  }
  if (recompute) {
    useFile <- FALSE
    if (is.null(sFile))
      sFile <- file.path(path, "L2Tprimees.rda")
  } else {
    useFile <- !is.null(sFile)
    if (useFile) {
      sFile <- file.path(path, sFile)
      if (!testPathForOutput(sFile, overwrite=FALSE))
        stop("the 'path' and/or 'sFile' provided are not correct or no permission to write in this file/path")  
        load(sFile) #-> 'mprime.approx'
    } else {
      MML2Tenv <- new.env()
      suppressWarnings(data("L2Tprimees", envir=MML2Tenv))
      if (exists("mprime.approx", envir=MML2Tenv, inherits=FALSE)) {
        mprime.approx <- with(MML2Tenv, mprime.approx)
        useFile <- TRUE
      } else {
        useFile <- FALSE
        if (is.null(sFile))
          sFile <- file.path(path, "L2Tprimees.rda")
        else
          sFile <- file.path(path, sFile)
      }
    }
  }
  
  if (!useFile) {
    mm.la <- mmprimeL2Tint(lambda)        
    mprime.approx <- splinefun(lambda, mm.la, method = "monoH.FC")
    e <- environment(mprime.approx)
    assign("lambda.max", max(lambda), envir=e)
    if (!testPathForOutput(sFile, overwrite=TRUE))
      stop("the 'path' and/or 'sFile' provided are not correct or no permission to write in this file/path")
    save(mprime.approx, file = sFile)
  }
  return(mprime.approx)
}

mmprimeL2T <- function(lam, mprime.approx) {
  la.max <- environment(mprime.approx)$lambda.max
  z <- ((m <- lam) <= la.max)
  m[z] <- mprime.approx(lam[z])
  if(any(i <- !z)) m[i] <- 1/(2*sqrt(lam[i]))
  m
}

wL2T <- function(mu, mprime.approx) {
  la.max <- environment(mprime.approx)$lambda.max
  ifelse(mu >= la.max, sqrt(mu)/2, mmprimeL2T(mu, mprime.approx)*mu)
}

if (FALSE) {
  library(numDeriv)
  m.approx <- mk.m_L2(lambda=seq(0,200,0.0005), sFile="L2TesPrecise.rda", recompute=TRUE)
  mprime.approx <- mk.mprime_L2(lambda=seq(0,200,0.0005), sFile="L2TprimeesPrecise.rda",  recompute=TRUE)
  lam <- seq(0.01,100,0.01)
  lamg <- grad(function(x) mmL2T(x, m.approx), x=lam)
  lamgg <- grad(function(x) mmL2Tint(x, lam.max=200), x=lam)
  lamm <- mmprimeL2T(lam, mprime.approx)
  plot(lam, lamm, type="l")
  lines(lam, lamg, col=2)
  lines(lam, lamgg, col=3)
  
  plot(lam, lamm-lamg, type="l")
  plot(lam, lamm-lamgg, type="l")
}
