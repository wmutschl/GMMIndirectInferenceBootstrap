# Plots for the slides
library(MASS)
library(rgl)
setwd("c:/temp")

# Bivariate density and contour plot
x1 <- mvrnorm(1000,c(0,0),matrix(c(2,1,1,2),2,2))
x2 <- mvrnorm(1000,c(2,2),matrix(c(2,-1.6,-1.6,2),2,2))
x <- rbind(x1,x2)
a <- kde2d(x,n=50)
par(mfrow=c(1,2))
persp(a,xlab="X",ylab="Y",zlab="density",main="3D plot")
contour(a,xlab="X",ylab="Y",main="Contour plot")
savePlot("bvdensity.pdf","pdf")

# Plot of likelihood function for uniform distribution
Like <- function(theta,x) {
  n <- length(x)
  L0 <- ifelse(theta<max(x),0,(1/theta)^n)
  }

windows(8,5)
par(mar=c(5,4,1,1))
x <- runif(20,0,4)
g <- seq(0,6,length=500)
plot(g,Like(g,x),xlab=expression(theta),ylab="likelihood",t="l")
savePlot("likeliuniform.pdf","pdf")

# Inconsistent estimation and forecasting
windows(8,5)
n <- 100
x <- runif(n,min=10,max=50)
u <- (x-mean(x))*2+rnorm(n,sd=5)
y <- 2+3*x+u
plot(x,y,main="Positive correlation between u and X")
abline(2,3)
text(locator(1),"true regression line")
