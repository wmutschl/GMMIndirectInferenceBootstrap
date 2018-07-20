# Using nonlinear transformations of instruments as
# additional instruments

library(MASS)
library(AER)

# Generate data: both exogenous variables are correlated with u
n <- 10000

Omega <- matrix(c(
2,0.3,0.5,0.7,
0.3,1,0.5,0.7,
0.5,0.5,1,0,
0.7,0.7,0,1),nrow=4,ncol=4)

R <- 100
Z <- matrix(NA,nrow=R,ncol=3)
ZZ <- matrix(NA,nrow=R,ncol=3)
for(r in 1:R) {
  dat <- mvrnorm(n,c(5,5,0,5),Omega)
  x1 <- dat[,1]
  x2 <- dat[,2]
  u <- dat[,3]
  w <- dat[,4] # instrument
  y <- 1+2*x1+3*x2+u

  #OLS is inconsistent
  summary(lm(y~x1+x2))
  ols <- lm(y~x1+x2)
  
  # IV from AER-package
  w2 <- w^2
  w3 <- w^3
  obj <- ivreg(y~x1+x2|w2+w3)
  
  #store estimates
  Z[r,] <- coefficients(obj)
  ZZ[r,] <- coefficients(ols)
  }

print(apply(Z,2,median))
print(apply(Z,2,mean))
print(apply(Z,2,sd))

print(apply(ZZ,2,median))
print(apply(ZZ,2,mean))
print(apply(ZZ,2,sd))

par(mfrow=c(3,1))
truehist(Z[,1])
truehist(Z[,2])
truehist(Z[,3])
