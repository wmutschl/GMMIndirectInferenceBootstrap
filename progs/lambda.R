# Verteilung des Momenten-Schätzers für die Exponentialverteilung

library(MASS)

# TEIL I: Ziehe R=1000 Stichproben der Länge n=8
#         und berechne jedesmal den Schätzer lambdahat

lambda <- 2
R <- 1000
v <- rep(NA,R)

for(i in 1:R)
  {
  x <- rexp(n=8,rate=lambda)
  lambdahat <- 1/mean(x)
  v[i] <- lambdahat
  hist(v,breaks=seq(0,12,length=30),ylim=c(0,250),main="Verteilung von lambdahat",col="purple")
  # das Histogramm wird bei jedem Schleifendurchlauf neu erzeugt, damit man die Entwicklung sieht
  abline(v=2,col="light grey",lwd=3)
  if(i<=10) # in den ersten 10 Durchläufen wird die Stichprobe angezeigt
    {
    text(9,seq(250,166,length=8),paste("x[",1:8,"] =",round(x,5)),pos=4)
    lines(c(9,12.5),c(158,158))
    text(9,146,paste("xbar =",round(mean(x),5)),pos=4)
    text(9,132,paste("lmbd=",round(1/mean(x),5)),pos=4)
    arrows(9,125,v[i]+0.25,10)
    locator(1) # das Programm wartet auf einen Mausklick (im Grafikfenster)
    }
  }
locator(1) # Mausklick
abline(v=mean(v),lwd=4) # Ergänzen des Durchschnitts (Schätzer für den Erwartungswert des Schätzers)
text(mean(v),230,paste("E(lambdahat)",round(mean(v),3)),pos=4)


# TEIL II (im wesentlichen wie in zgs.R)

mu <- 1/lambda
sigma <- sqrt(1/lambda^2)

R <- 10000
g <- seq(-2.5,2.5,length=500)

for(n in 1:50)
  {
  x <- matrix(rexp(n*R,rate=lambda),R,n)
  xbar <- apply(x,1,mean)
  z <- sqrt(n)*(xbar-mu)
  truehist(z,xlim=c(-2.5,2.5), ylim=c(0,2), main=paste("n = ",n))
  lines(g,dnorm(g,mean=0,sd=sigma),col="red",lwd=3)
  if(n<3) locator(1)
  }
