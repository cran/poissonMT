#############################################################
#	functions for poissonMT
#	Author: C. Agostinelli, M. Valdora and V.J. Yohai
#	Maintainer e-mail: claudio.agostinelli@unitn.it
#	Date: 17 October 2018
#	Version: 0.2-4
#	Copyright (C) 2018 C. Agostinelli M. Valdora and V.J. Yohai
#############################################################

### Taken from robustbase version 0.92-4
mk.m_rho <- function(cw, opt.method = c("L-BFGS-B", "Brent",
    "Nelder-Mead", "BFGS", "CG", "SANN"), lambda = c(seq(0,2.9,
    by=0.1), seq(3,100)), reltol = sqrt(.Machine$double.eps),
    trace = 0, sFile = NULL,
    path = poissonMTgetwd(),  
    recompute  = getOption("poissonMT:m_rho_recompute", FALSE)) {
  if (is.null(path))
    stop("please set up the working directory for this package by calling 'poissonMTsetwd'")  
  if (recompute) {
    useFile <- FALSE
    if (is.null(sFile))
      sFile <- file.path(path, paste0("MTes_", format(cw), ".rda"))
  } else {
    useFile <- !is.null(sFile)
    if (useFile) {
      sFile <- file.path(path, sFile)
      if (!testPathForOutput(sFile, overwrite=FALSE))
        stop("the 'path' and/or 'sFile' provided are not correct or no permission to write in this file/path")
        load(sFile) #-> 'm.approx
    ## check if its cw was very close to this one
      if (cw.ok <- is.numeric(cw0 <- environment(m.approx)$cw))
        cw.ok <- (abs(cw - cw0) < 0.001)  
    } else {
      MMMTenv <- new.env()
      suppressWarnings(data(list=paste0("MTes_", format(cw)), envir=MMMTenv))
      if (exists("m.approx", envir=MMMTenv, inherits=FALSE)) {
        m.approx <- with(MMMTenv, m.approx)
        useFile <- TRUE
    ## check if its cw was very close to this one
        if (cw.ok <- is.numeric(cw0 <- environment(m.approx)$cw))
          cw.ok <- (abs(cw - cw0) < 0.001)
      } else {
        useFile <- FALSE
        if (is.null(sFile))
          sFile <- file.path(path, paste0("MTes_", format(cw), ".rda"))
        else
          sFile <- file.path(path, sFile)
      }
    }
  }

  if (!useFile || !cw.ok) {
    nl <- length(lambda)
    mm.la <- numeric(nl)
    s.la <- sqrt(lambda)

    ## MM: Speedwise,   Brent > L-BFGS-B > BFGS > ..  for cw >= ~ 1.5
    ##         L-BFGS-B > Brent  for  cw ~= 1
    opt.method <- match.arg(opt.method)
    oCtrl <- list(reltol=reltol, trace=trace)
    if (opt.method %in% c("Brent", "L-BFGS-B")) { ## use bounds
      if (opt.method == "L-BFGS-B")# yuck! why is this necessary!!
	oCtrl <- list(factr = 1/(10*reltol), trace=trace)
      for (i in seq_len(nl))
        mm.la[i] <- optim(s.la[i], espRho, lam = lambda[i], cw = cw,
         method = opt.method, control = oCtrl,
         lower = 0, upper = .01 + 2*s.la[i])$par
    } else {
      for (i in seq_len(nl))
        mm.la[i] <- optim(s.la[i], espRho, lam = lambda[i], cw = cw,
	method = opt.method, control = oCtrl)$par
    }
    m.approx <- splinefun(lambda, mm.la, method = "monoH.FC")
    e <- environment(m.approx)
    assign("lambda.max", max(lambda), envir=e)
    assign("cw", cw, envir=e)
    if (!testPathForOutput(sFile, overwrite=TRUE))
      stop("the 'path' and/or 'sFile' provided are not correct or no permission to write in this file/path")
    save(m.approx, file = sFile)
  }
  m.approx
}

mm <- function (lam, m.approx) {
  la.max <- environment(m.approx)$lambda.max
  z <- ((m <- lam) <= la.max)
  m[z] <- m.approx(lam[z])
  if (any(i <- !z)) 
    m[i] <- sqrt(lam[i])
  m
}

espRho <- function (lam, xx, cw) {
  k <- seq(as.integer((max(0, xx - cw))^2), as.integer((xx + cw)^2) + 1L)
  inner <- (rhoS.k <- Mchi(x=sqrt(k) - xx, cc=cw, psi="bisquare")) < 1
  ii <- k[inner]
  terminos <- rhoS.k[inner] * dpois(ii, lam)
  if ((len.ii <- length(ii)) > 0) {
    primero <- ii[1]
    ultimo <- ii[len.ii]
    ppois(primero - 1, lam) + sum(terminos) + ppois(ultimo, lam, lower.tail = FALSE)
  } else 1
}
  
mmd.int <- function(lam,cw,alpha,m.approx) {
  qq1 <- qpois(alpha,lam)
  qq2 <- qpois(1-alpha,lam)
  ind <- qq1:qq2
  k.. <- sqrt(ind) - mm(lam, m.approx)
  dP <- dpois(ind,lam)
  rr1 <- (-dP+(ind*dP/lam)) * Mpsi(x=k.., cc=cw, psi="bisquare", deriv=0)
  rr2 <- dP*Mpsi(x=k.., cc=cw, psi="bisquare", deriv=1)
  rr1 <- sum(rr1)
  rr2 <- sum(rr2)
  rr1/rr2
}

mmdd.int <- function (lam, cw, alpha, m.approx) {
  qq1 <- qpois(alpha,lam)
  qq2 <- qpois(1-alpha,lam)
  ind <- qq1:qq2
  k.. <- sqrt(ind) - mm(lam, m.approx)
  dP <- dpois(ind,lam)
  NUM <- (-dP+(ind*dP/lam)) * Mpsi(x=k.., cc=cw, psi="bisquare", deriv=0)
  DEN <- dP * Mpsi(x=k.., cc=cw, psi="bisquare", deriv=1)
  NUM <- sum(NUM)
  DEN <- sum(DEN)
  mm1 <- NUM/DEN
  NUMP <- ddpois(ind, lam) *
    Mpsi(x=k.., cc=cw, psi="bisquare", deriv=0) - (-dP + (ind * dP/lam)) * 
    Mpsi(x=k.., cc=cw, psi="bisquare", deriv=1) * mm1
  DENP <- (-dP + (ind * dP/lam)) *
    Mpsi(x=k.., cc=cw, psi="bisquare", deriv=1) - dP *
    Mpsi(x=k.., cc=cw, psi="bisquare", deriv=2) * mm1
  NUMP <- sum(NUMP)
  DENP <- sum(DENP)
  list(mm1=mm1, mm2=(NUMP * DEN - DENP * NUM)/DEN^2)
}

ddpois <- function(x, lam) {
    ## The second derivative of the Poisson probability function
    dpois(x,lam)*(1-(2*x/lam)+((x^2)/(lam^2))-(x/(lam^2)))
}
  
mk.mprime_rho <- function(lambda = c(seq(0,2.9, by=0.1), seq(3,200)),
    alpha = 0.00001, cw=2.3, sFile = NULL,
    path = poissonMTgetwd(),  
    recompute  = getOption("poissonMT:m_rho_recompute", FALSE)) {
  if (is.null(path)) {
    stop("please set up the working directory for this package by calling 'poissonMTsetwd'")
  }
  if (recompute) {
    useFile <- FALSE
    if (is.null(sFile))
      sFile <- file.path(path, paste0("MTesPrime_", format(cw), "-", format(alpha), ".rda"))
  } else {
    useFile <- !is.null(sFile)
    if (useFile) {
      sFile <- file.path(path, sFile)
      if (!testPathForOutput(sFile, overwrite=FALSE))
        stop("the 'path' and/or 'sFile' provided are not correct or no permission to write in this file/path")
      load(sFile) #-> 'mprime.approx
    ## check if its cw was very close to this one
      if (cw.ok <- is.numeric(cw0 <- environment(mprime.approx)$cw))
        cw.ok <- (abs(cw - cw0) < 0.001)  
    ## check if its alpha was very close to this one
        if (alpha.ok <- is.numeric(alpha0 <- environment(mprime.approx)$alpha))
          alpha.ok <- (abs(alpha - alpha0)/alpha0 < 0.001)
    } else {
      MMMTprimeenv <- new.env()
      suppressWarnings(data(list=paste0("MTesPrime_", format(cw), "-", format(alpha)), envir=MMMTprimeenv))
      if (exists("mprime.approx", envir=MMMTprimeenv, inherits=FALSE)) {
        mprime.approx <- with(MMMTprimeenv, mprime.approx)
        useFile <- TRUE
    ## check if its cw was very close to this one
        if (cw.ok <- is.numeric(cw0 <- environment(mprime.approx)$cw))
          cw.ok <- (abs(cw - cw0) < 0.001)
    ## check if its alpha was very close to this one
        if (alpha.ok <- is.numeric(alpha0 <- environment(mprime.approx)$alpha))
          alpha.ok <- (abs(alpha - alpha0)/alpha0 < 0.001)
      } else {
        useFile <- FALSE
        if (is.null(sFile))
          sFile <- file.path(path, paste0("MTesPrime_", format(cw), "-", format(alpha), ".rda"))
        else
          sFile <- file.path(path, sFile)
      }
    }
  }

  if (!useFile || !cw.ok || !alpha.ok) {    
    m.approx <- mk.m_rho(cw=cw)  
    mmd.la <- sapply(lambda, function(x) mmd.int(x,cw,alpha,m.approx))
    mprime.approx <- splinefun(lambda, mmd.la, method = "monoH.FC")
    e <- environment(mprime.approx)
    assign("lambda.max", max(lambda), envir=e)
    assign("cw", cw, envir=e)
    assign("alpha", alpha, envir=e)
    if (!testPathForOutput(sFile, overwrite=TRUE))
      stop("the 'path' and/or 'sFile' provided are not correct or no permission to write in this file/path")
    save(mprime.approx, file = sFile)
  }
  mprime.approx
}

mmd <- function(lam, mprime.approx) {
  la.max <- environment(mprime.approx)$lambda.max
  z <- ((m <- lam) <= la.max)
  m[z] <- mprime.approx(lam[z])
  if(any(i <- !z)) m[i] <- 1/(2*sqrt(lam[i]))
  m
}
