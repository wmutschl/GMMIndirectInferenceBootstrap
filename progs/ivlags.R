# IV estimation with lagged variables

library(AER)
set.seed(123)
TT <- 50000

# Example 1: Measurement errors in time series

xstar <- as.numeric(filter(rnorm(TT,sd=sqrt(1-0.6^2)),0.6,method="r",init=rnorm(1)))
v <- rnorm(TT,sd=0.5)
x <- xstar+v
u <- rnorm(TT)
# True values
alpha <- 2
beta <- 3
y <- alpha+beta*xstar+u
# OLS is inconsistent
print(lm(y~x))
# IV is consistent
print(ivreg(y[-1]~x[-1]|x[-TT]))


# Example 2: Omitted variables in time series

x1 <- rep(NA,TT)
x2 <- rep(NA,TT)
x1[1] <- rnorm(1)
x2[1] <- rnorm(1)
for(tt in 2:TT) {
  x1[tt] <- 0.6*x1[tt-1]+0.3*x2[tt-1]+rnorm(1,sd=0.3)
  x2[tt] <- 0.3*x1[tt-1]+0.6*x2[tt-1]+rnorm(1,sd=0.3)
}
# True values: 
alpha <- 1
beta1 <- 2
beta2 <- 3
y <- alpha+beta1*x1+beta2*x2+rnorm(TT)
# OLS is inconsistent
print(lm(y~x1))
# IV is also inconsistent
print(ivreg(y[-1]~x1[-1]|x1[-TT]))


# Example 3: Endogeneity in time series

x <- rep(NA,TT)
y <- rep(NA,TT)
x[1] <- rnorm(1,mean=10)
y[1] <- rnorm(1,mean=10)
u <- rnorm(TT,sd=0.3)
v <- rnorm(tt,sd=0.3)
# True values: 
alpha <- 1
beta1 <- 0.6
beta2 <- 0.3 
gamma <- 1
delta1 <- 0.6
delta2 <- 0.3
for(tt in 2:TT) {
  c0 <- 1-beta1*delta1
  y[tt] <- (alpha+beta1*gamma+beta1*delta2*x[tt-1]+beta2*y[tt-1]+beta1*v[tt]+u[tt])/c0
  x[tt] <- gamma+delta1*y[tt]+delta2*x[tt-1]+v[tt]
}
# OLS is inconsistent
print(lm(y[-1]~x[-1]+y[-TT]))
# IV is consistent
print(ivreg(y[-1]~x[-1]+y[-TT]|x[-TT]+y[-TT]))
