# Illustration of errors in variables

library(MASS)

TT <- 500
beta0 <- c(1,2,3,-4)
Xstar <- cbind(1,mvrnorm(TT,rep(5,3),matrix(0.6,3,3)+0.4*diag(3)))
y <- Xstar%*%beta0+rnorm(TT,sd=0.5)
# Regression on true matrix Xstar
lm(y~Xstar-1)
solve(t(Xstar)%*%Xstar)%*%t(Xstar)%*%y

X <- Xstar
X[,4] <- Xstar[,4]+rnorm(TT,sd=.5)
# Regression on measured matrix X
lm(y~X-1)
solve(t(X)%*%X)%*%t(X)%*%y
