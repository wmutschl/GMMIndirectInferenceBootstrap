# Finite sample properties of IV estimators
library(AER)
library(MASS)
library(evir)
set.seed(123)

# Example 1: Measurement errors in time series

# True values
alpha <- 2
beta <- 3

TT <- 50
R <- 10000
Z <- rep(NA,R)
for(r in 1:R) {
  xstar <- as.numeric(filter(rnorm(TT,sd=sqrt(1-0.6^2)),0.6,method="r",init=rnorm(1)))
  x <- xstar+rnorm(TT,sd=0.5)
  u <- rnorm(TT)
  y <- alpha+beta*xstar+u
  w <- x[-TT]
  #obj <- ivreg(y[-1]~x[-1]|w)
  #betahat <- obj$coefficients[2]
  betahat <- sum((w-mean(w))*(y[-1]-mean(y[-1])))/sum((w-mean(w))*(x[-1]-mean(x[-1])))
  Z[r] <- betahat
}

# Histogram over entire range
truehist(Z,main="Distribution of IV estimator")

# Descriptive statistics
print(paste("Range   =",range(Z)))
print(paste("Mean    =",mean(Z)))
print(paste("Std.dev.=",sd(Z)))
print(paste("Median  =",median(Z)))

# Main part of the histogram plus density of fitted normal
truehist(Z[abs(Z-3)<10],main="Distribution of IV estimator")
g <- seq(-7,13,length=500)
lines(g,dnorm(g,mean(Z),sd(Z)),lwd=2)

# Hill plot indicates that the tail index is unity (ie. no first moment)
hill(Z,end=40)
abline(h=1,lty="dashed")


# Exercise 13.5
# (Davidson and MacKinnon, 2004, exercise 8.10 and 8.11)

TT <- 10
beta0 <- 0
pi0 <- 2
sigma.u <- sigma.v <- 1
rho <- 0.5
R <- 1000
Z <- matrix(NA,R,2)
for(r in 1:R) {
  w <- rnorm(TT)
  w <- w/sqrt(w%*%w/TT)
  uv <- mvrnorm(TT,c(0,0),matrix(c(1,rho,rho,1),2,2))
  x <- pi0*w+sigma.v*uv[,2]
  y <- beta0*x+sigma.u*uv[,1]
  Z[r,1] <- coefficients(lm(y~x-1))
  Z[r,2] <- coefficients(ivreg(y~x-1|w))
}
# OLS estimator
plot(ecdf(Z[,1]),main="Distributions of estimators")
# Simple IV estimator
lines(ecdf(Z[,2]),col="red")
# Asymptotic distribution
g <- seq(-1,1.5,length=500)
lines(g,pnorm(g,sd=sigma.u/sqrt(TT*pi0^2)),col="green")
