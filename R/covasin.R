### Taken from robustbase 0.93.1 and modified
covasin <- function(x,y,beta,cw, m.approx,w)
{
    p <- ncol(x)
    n <- length(y)
    mm1 <- mm2 <- numeric(n)
    de <- nu <- matrix(0,p,p)

    lam <- x%*%beta
    elam <- exp(lam)
    r <- sqrt(y) - mm(elam, m.approx)
    psi0 <- Mpsi(x=r, cc=cw, psi="bisquare", deriv=0)
    psi1 <- Mpsi(x=r, cc=cw, psi="bisquare", deriv=1)
    for ( i in 1:n) {
      temp <- mmdd.int(elam[i], cw, 0.001, m.approx)
      mm1[i] <- temp$mm1
      mm2[i] <- temp$mm2
    }
    nu1 <- w*psi0*mm1*elam
    de1 <- -psi1*(mm1^2)*(elam^2)+psi0*mm2*(elam^2)+psi0*mm1*elam
    de1 <- w*de1

    for (i in 1:n) { ## FIXME (?)  -- can be vectorized
        zzt <- tcrossprod(x[i,])
        nu <- nu+ (nu1[i]^2)*zzt
        de <- de+ de1[i]*zzt
    }
    nu <- nu/n
    de <- solve(de/n)
    ## Cov_{asympt.} =
    de %*% nu %*% t(de) / n
}
