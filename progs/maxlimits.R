library(MASS)

# Limits of maxima (I)
n <- 500
N <- 10000
M <- rep(NA,N)
for(i in 1:N) {
  x <- rnorm(n)
  M[i] <- max(x)
  }

dn <- sqrt(2*log(n))-(log(4*pi)+log(log(n)))/(2*sqrt(2*log(n)))
cn <- (2*log(n))^(-0.5)
R <- (M-dn)/cn
g <- seq(min(R),max(R),length=400)
truehist(R)
lines(g,exp(-g-exp(-g)))


# Limits of maxima (II)
n <- 500
N <- 10000
M <- rep(NA,N)
for(i in 1:N) {
  x <- rt(n,df=1.5)
  M[i] <- max(x)
  }

cn <- qt(1-1/n,df=1.5)
R <- M/cn
g <- seq(min(R),max(R),length=400)
truehist(R)
lines(g,1.5*g^(-2.5)*exp(-g^(-1.5)))


# Limits of maxima (III)
n <- 500
N <- 10000
M <- rep(NA,N)
for(i in 1:N) {
  x <- runif(n)
  M[i] <- max(x)
  }

dn <- 1
cn <- 1/n
R <- (M-dn)/cn
g <- seq(min(R),max(R),length=400)
truehist(R)
lines(g,exp(g))
