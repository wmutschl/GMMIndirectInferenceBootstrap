library(MASS)

g <- seq(-5,5,length=300)

for(n in 1:9)
  {
  x <- matrix(rcauchy(10000*n),10000,n)
  xbar <- apply(x,MARGIN=1,mean)
  plot(density(xbar,from=-5,to=5),ylim=c(0,0.4),main=paste("n=",n))
  lines(g,dcauchy(g),lwd=3,col="pink")
  }
  
