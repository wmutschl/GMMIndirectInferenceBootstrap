# Illustration eines stochastischen Prozesses

wait <- function(hsec)
  {
  q0 <- proc.time()[3]
  while((proc.time()[3]-q0)<hsec/100) z <- 0
  }
  
p0 <- proc.time()[3]
while(proc.time()[3]-p0<10)
  {
  plot(10*exp(cumsum(rnorm(250,sd=0.25/sqrt(250)))),type="l",ylim=c(0,20),xlab="Zeit t",ylab="X(t)",main="Zeit variabel, omega variabel")
  wait(8)
  }

plot(10*exp(cumsum(rnorm(250,sd=0.25/sqrt(250)))),type="l",ylim=c(0,20),xlab="Zeit t",ylab="X(t)",main="Zeit variabel, omega fest")
locator(1)

p0 <- proc.time()[3]
while(proc.time()[3]-p0<10)
  {
  x <- 10*exp(cumsum(rnorm(250,sd=0.25/sqrt(250))))
  plot(x,type="l",ylim=c(0,20),xlab="Zeit t",ylab="X(t)",main="Zeit fest, omega variabel")
  abline(v=150,lty="dashed")
  points(150,x[150],col="red",lwd=5,pch=19)
  wait(8)
  }

x <- 10*exp(cumsum(rnorm(250,sd=0.25/sqrt(250))))
plot(x,type="l",ylim=c(0,20),xlab="Zeit t",ylab="X(t)",main="Zeit fest, omega fest")
abline(v=150,lty="dashed")
points(150,x[150],col="red",lwd=5,pch=19)
