## from robustbase version 0.93.1
robXweights <- function (wts, X, intercept = TRUE) {
    stopifnot(length(d <- dim(X)) == 2, is.logical(intercept))
    nobs <- d[1]
    if (d[2]) {
        if (is.character(wts)) {
            switch(wts, none = rep.int(1, nobs), hat = wts_HiiDist(X)^2, 
                robCov = wts_RobDist(X, intercept, covFun = cov.rob), 
                covMcd = wts_RobDist(X, intercept, covFun = covMcd), 
                stop("Weighting method", sQuote(wts), " is not implemented"))
        }
        else if (is.list(wts)) {
            if (length(wts) == 1 && is.function(covF <- wts[[1]])) 
                wts_RobDist(X, intercept, covFun = covF)
            else stop("if a list, weights.on.x must contain a covariance function such as covMcd()")
        }
        else if (is.function(wts)) {
            wts(X, intercept)
        }
        else {
            if (!is.numeric(wts) || length(wts) != nobs) 
                stop(gettextf("weights.on.x needs %d none-negative values", 
                  nobs), domain = NA)
            if (any(wts) < 0) 
                stop("All weights.on.x must be none negative")
        }
    }
    else rep.int(1, nobs)
}

wts_HiiDist <- function(X) {
  ## Hii := diag( tcrossprod( qr.Q(qr(X)) ) ) == rowSums( qr.Q(qr(X)) ^2 ) :
  x <- qr(X)
  Hii <- rowSums(qr.qy(x, diag(1, nrow = NROW(X), ncol = x$rank))^2)
  (1-Hii)
}

wts_RobDist <- function(X, intercept, covFun) {
  D2 <- if (intercept) { ## X[,] has intercept column which should not be used for rob.wts
    X <- X[, -1, drop=FALSE]
    Xrc <- covFun(X)
    mahalanobis(X, center = Xrc$center, cov = Xrc$cov)
  } else { ## X[,]  can be used directly
    if (!is.matrix(X)) X <- as.matrix(X)
    Xrc <- covFun(X)
    S <- Xrc$cov + tcrossprod(Xrc$center)
    mahalanobis(X, center = FALSE, cov = S)
  }
  p <- ncol(X) ## E[chi^2_p] = p
  1/sqrt(1+ pmax.int(0, 8*(D2 - p)/sqrt(2*p)))
}

    