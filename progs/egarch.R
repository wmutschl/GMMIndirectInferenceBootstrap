# Die Verteilung des Autokorrelationskoeffizienten in einem t-GARCH(1,1)-Prozess

library(MASS)
T <- 2500
alpha0 <- 0.05
alpha1 <- 0.1
beta1 <- 0.8
nu <- 7
sgm2 <- 1


# Funktion zum Generieren eines Pfads des t-GARCH(1,1)-Prozesses
garch.sim <- function(T,alpha0,alpha1,beta1,nu,sgm2)
  {
  e <- rt(T+100,df=nu)*sqrt((nu-2)/nu)
  x <- rep(0,T+100)
  for(i in 2:(T+100))
    {
    x[i] <- sqrt(sgm2)*e[i]
    sgm2 <- alpha0+alpha1*x[i-1]^2+beta1*sgm2
    }
  return(x[101:(T+100)])
  }

v <- rep(NA,1000)
for(i in 1:1000)
  {
  x <- tgarch.sim(T,alpha0,alpha1,beta1,nu,sgm2)
  v[i] <- cor(x[-1],x[-T])
  }

truehist(v)
