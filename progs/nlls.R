# Nonlinear least squares estimation

# Exercise 11.1.1

# Load the dataset
dat <- read.csv("expgrowth.csv")

# Definition of the objective function
# The first argument must be the parameter vector
squarediffs <- function(param,dat) {
  alpha <- param[1]
  beta <- param[2]
  d <- dat$y-exp(alpha+beta*dat$x)
  return(sum(d^2))
}

# Minimize the objective function
obj <- optim(c(1,0),squarediffs,dat=dat)
estimates <- obj$par
alphahat <- estimates[1]
betahat <- estimates[2]

# Plot the data and the estimated regression curve
plot(dat$x,dat$y)
g <- seq(0,40,length=500)
lines(g,exp(alphahat+betahat*g))


# Exercise 11.1.2

# Load the dataset
dat <- read.csv("DMacK1.csv")

# Definition of the objective function
# The first argument must be the parameter vector
squarediffs <- function(param,dat) {
  beta1 <- param[1]
  beta2 <- param[2]
  d <- dat$y-(beta1+beta2*dat$x1+1/beta2*dat$x2)
  return(sum(d^2))
}

# Minimize the objective function
obj <- optim(c(1,1),squarediffs,dat=dat)
estimates <- obj$par
beta1hat <- estimates[1]
beta2hat <- estimates[2]

# Plot the data and the estimated regression plane
library(rgl)
plot3d(dat$x1,dat$x2,dat$y)
gx1 <- seq(min(dat$x1),max(dat$x1),length=60)
gx2 <- seq(min(dat$x2),max(dat$x2),length=60)
yhat <- outer(gx1,gx2,function(x1,x2) beta1hat+beta2hat*x1+1/beta2hat*x2)
surface3d(gx1,gx2,yhat,col="light green")
