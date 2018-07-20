# Numerical estimation of the parameters of N(mu,sigma2)

# Generate sample
set.seed(1234)
x <- rnorm(n=50,mean=5,sd=3)

# Definition of loglikelihood
neglogl <- function(theta,x) {
  mu <- theta[1]
  sigma2 <- theta[2]
  return(-sum(log(dnorm(x,mean=mu,sd=sqrt(sigma2)))))
  }

# Optimization
obj <- optim(c(0,1),neglogl,x=x,hessian=TRUE)

# Point estimates
print(obj$par)

# Numerical covariance matrix
print(solve(obj$hessian))
