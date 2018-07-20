# Programm zur Illustration der Konsistenz eines Schätzers

R <- 500
g <- seq(-3,3,length=250)

for(n in c(1:50,seq(100,1000,by=50)))
  {
  x <- matrix(rnorm(n*R),R,n)
  xbar <- apply(x,1,mean)+2/n
  truehist(xbar,xlim=c(-3,3), ylim=c(0,15), main=paste("n = ",n))
  lines(g,dnorm(g-2/n,sd=sqrt(1/n)),col="red")
  if(n<6) locator(1) # Warten auf Mausklick
  }
