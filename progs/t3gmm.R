library(gmm)
library(MASS)
library(mvtnorm)

g <- function(theta,x) {
  beta <- theta[1]
  sigma <- theta[2]
  u <- x[,2]-beta*x[,1]
  gmat <- cbind(u,u*x[,1],u^2-sigma^2)
  return(gmat)
  }

R <- 1000
n <- 5000
nu <- 30
Z <- matrix(NA,R,2)
for(r in 1:R) {
  x1 <- rt(n,df=nu)/sqrt(nu/(nu-2))
  u <- rt(n,df=nu)/sqrt(nu/(nu-2))
  x2 <- 0.9*x1+u
  a <- gmm(g,x=cbind(x1,x2),t0=c(0.9,1),wmatrix="optimal")
  Z[r,] <- (abs(coefficients(a))-c(0.9,1))/sqrt(diag(vcov(a)))
  }
apply(Z,2,mean)

truehist(Z[,2])
g <- seq(min(Z[,2]),max(Z[,2]),length=500)
lines(g,dnorm(g,mean(Z[,2]),sd(Z[,2])))
