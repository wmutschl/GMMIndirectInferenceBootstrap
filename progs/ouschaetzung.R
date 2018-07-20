# Illustration zur Schätzung eines OU-Prozesses

sim.OU <- function(eta,sigma=1,meany=0,N=1000,T=1,y0=0) {
  dd <- T/N
  d <- seq(0,T,length=N)
  y <- filter(rnorm(N-1,mean=meany-exp(-eta*dd)*meany,sd=sqrt(sigma^2/(2*eta)*(1-exp(-2*dd*eta)))),exp(-eta*dd),method="recursive",init=y0)
  return(list(x=d,y=c(y0,y)))
  }

delay <- function(c) {
  p0 <- proc.time()[3]
  while(proc.time()[3]-p0<c) {}
  }
  
# Der diskret beobachtete Prozess
x <- sim.OU(2,sigma=1,meany=4,T=50,N=1001,y0=4)
plot(x,type="l",xlab="Zeit t",ylab="X(t)",main="Realisierter Prozess (Pfad)",ylim=c(0,10))
a <- locator(1)
plot(x,type="l",xlab="Zeit t",ylab="X(t)",main="Beobachteter Prozess (Pfad)",ylim=c(0,10))
b <- list(x=x$x[seq(1,1001,by=20)],y=x$y[seq(1,1001,by=20)])
lines(x,col="grey")
points(b,col="red",pch=19,lwd=2)
lines(b,col="red",lwd=2)
a <- locator(1)

# Passt ein Prozess mit anderem Mittelwert ?
for(i in 1:25) {
  plot(b,type="l",xlab="Zeit t",ylab="X(t)",main="Variation von alpha",ylim=c(0,10),col="red",lwd=2)
  text(25,10,"alpha=6, eta=2, sigma=1")
  points(b,col="red",pch=19,lwd=2)
  y <- sim.OU(2,sigma=1,meany=6,T=50,N=1001,y0=6)
  lines(y,col="grey")
  yb <- list(x=y$x[seq(1,1001,by=20)],y=y$y[seq(1,1001,by=20)])
  points(yb,pch=19,lwd=2)
  lines(yb,lwd=2)
  delay(0.1)
  }
a <- locator(1)
for(i in 1:25) {
  plot(b,type="l",xlab="Zeit t",ylab="X(t)",main="Variation von alpha",ylim=c(0,10),col="red",lwd=2)
  text(25,10,"alpha=3, eta=2, sigma=1")
  points(b,col="red",pch=19,lwd=2)
  y <- sim.OU(eta=2,sigma=1,meany=3,T=50,N=1001,y0=3)
  lines(y,col="grey")
  yb <- list(x=y$x[seq(1,1001,by=20)],y=y$y[seq(1,1001,by=20)])
  points(yb,pch=19,lwd=2)
  lines(yb,lwd=2)
  delay(0.1)
  }
a <- locator(1)
for(i in 1:25) {
  plot(b,type="l",xlab="Zeit t",ylab="X(t)",main="Variation von alpha",ylim=c(0,10),col="red",lwd=2)
  text(25,10,"alpha=4, eta=2, sigma=1")
  points(b,col="red",pch=19,lwd=2)
  y <- sim.OU(eta=2,sigma=1,meany=4,T=50,N=1001,y0=4)
  lines(y,col="grey")
  yb <- list(x=y$x[seq(1,1001,by=20)],y=y$y[seq(1,1001,by=20)])
  points(yb,pch=19,lwd=2)
  lines(yb,lwd=2)
  delay(0.1)
  }
a <- locator(1)

# Passt ein Prozess mit einem anderen sigma ?
for(i in 1:25) {
  plot(b,type="l",xlab="Zeit t",ylab="X(t)",main="Variation von sigma",ylim=c(0,10),col="red",lwd=2)
  text(25,10,"alpha=4, eta=2, sigma=3")
  points(b,col="red",pch=19,lwd=2)
  y <- sim.OU(eta=2,sigma=3,meany=4,T=50,N=1001,y0=4)
  lines(y,col="grey")
  yb <- list(x=y$x[seq(1,1001,by=20)],y=y$y[seq(1,1001,by=20)])
  points(yb,pch=19,lwd=2)
  lines(yb,lwd=2)
  delay(0.1)
  }
a <- locator(1)
for(i in 1:25) {
  plot(b,type="l",xlab="Zeit t",ylab="X(t)",main="Variation von sigma",ylim=c(0,10),col="red",lwd=2)
  text(25,10,"alpha=4, eta=2, sigma=0.1")
  points(b,col="red",pch=19,lwd=2)
  y <- sim.OU(eta=2,sigma=0.1,meany=4,T=50,N=1001,y0=4)
  lines(y,col="grey")
  yb <- list(x=y$x[seq(1,1001,by=20)],y=y$y[seq(1,1001,by=20)])
  points(yb,pch=19,lwd=2)
  lines(yb,lwd=2)
  delay(0.1)
  }
a <- locator(1)
for(i in 1:25) {
  plot(b,type="l",xlab="Zeit t",ylab="X(t)",main="Variation von sigma",ylim=c(0,10),col="red",lwd=2)
  text(25,10,"alpha=4, eta=2, sigma=1")
  points(b,col="red",pch=19,lwd=2)
  y <- sim.OU(eta=2,sigma=1,meany=4,T=50,N=1001,y0=4)
  lines(y,col="grey")
  yb <- list(x=y$x[seq(1,1001,by=20)],y=y$y[seq(1,1001,by=20)])
  points(yb,pch=19,lwd=2)
  lines(yb,lwd=2)
  delay(0.1)
  }
a <- locator(1)

# Passt ein Prozess mit einem anderen eta ?
for(i in 1:25) {
  plot(b,type="l",xlab="Zeit t",ylab="X(t)",main="Variation von eta (und sigma)",ylim=c(0,10),col="red",lwd=2)
  text(25,10,"alpha=4, eta=20, sigma=3.162")
  points(b,col="red",pch=19,lwd=2)
  y <- sim.OU(eta=20,sigma=3.162,meany=4,T=50,N=1001,y0=4)
  lines(y,col="grey")
  yb <- list(x=y$x[seq(1,1001,by=20)],y=y$y[seq(1,1001,by=20)])
  points(yb,pch=19,lwd=2)
  lines(yb,lwd=2)
  delay(0.1)
  }
a <- locator(1)
for(i in 1:25) {
  plot(b,type="l",xlab="Zeit t",ylab="X(t)",main="Variation von eta (und sigma)",ylim=c(0,10),col="red",lwd=2)
  text(25,10,"alpha=4, eta=0.1, sigma=0.224")
  points(b,col="red",pch=19,lwd=2)
  y <- sim.OU(eta=0.1,sigma=0.224,meany=4,T=50,N=1001,y0=4)
  lines(y,col="grey")
  yb <- list(x=y$x[seq(1,1001,by=20)],y=y$y[seq(1,1001,by=20)])
  points(yb,pch=19,lwd=2)
  lines(yb,lwd=2)
  delay(0.1)
  }
a <- locator(1)
for(i in 1:25) {
  plot(b,type="l",xlab="Zeit t",ylab="X(t)",main="Variation von eta (und sigma)",ylim=c(0,10),col="red",lwd=2)
  text(25,10,"alpha=4, eta=2, sigma=1")
  points(b,col="red",pch=19,lwd=2)
  y <- sim.OU(eta=2,sigma=1,meany=4,T=50,N=1001,y0=4)
  lines(y,col="grey")
  yb <- list(x=y$x[seq(1,1001,by=20)],y=y$y[seq(1,1001,by=20)])
  points(yb,pch=19,lwd=2)
  lines(yb,lwd=2)
  delay(0.1)
  }
