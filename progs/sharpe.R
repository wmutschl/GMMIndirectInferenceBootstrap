# Zahlenbeispiel für die asymptotische Verteilung der Sharpe-Ratio

library(MASS)

# TEIL I: Ziehen einer Stichprobe, Berechnung des 0.95-Konfidenzintervalls

# Die Tages-Renditen seien iid N(0.034,2.5)
# Der risikolose Tages-Zinssatz sei 0.014

mu <- 0.034
sigma2 <- 2.5
Rf <- 0.014

# Berechnung der wahren SR
SR <- (mu-Rf)/sqrt(sigma2)

# Generieren eines Renditevektors der Länge n=1000
n <- 1000
x <- rnorm(n,mean=mu,sd=sqrt(sigma2))

# Berechnung des Schätzwerts für SR
muhat <- mean(x)
sigmahat <- sd(x)
SRhat <- (muhat-Rf)/sigmahat

# Berechnen des Konfidenzintervalls nach der Delta-Methode für die Stichprobe x

# 1. Berechnung des Gradienten (an der Stelle der Schätzwerte)
DSR <- matrix(c(1/sigmahat,-0.5*(muhat-Rf)/(sigmahat^3)),1,2)

# 2. Berechnung der (geschätzten) Kovarianzmatrix
Sigma <- matrix(c(sigmahat^2,0,0,2*sigmahat^4),2,2)

# 3. Berechnung der (geschätzten) Varianz
VSR <- DSR%*%Sigma%*%t(DSR)

# 4. Berechnung des 0.95-Konfidenzintervalls
print(SRhat-1.96*sqrt(VSR/n))
print(SRhat+1.96*sqrt(VSR/n))



# TEIL II: Die Verteilung von SRhat

# SRhat ist eine Zufallsvariable ! 
# Simulation ihrer Verteilung

v <- rep(NA,10000)
for(i in 1:10000)
  {
  x <- rnorm(n,mean=mu,sd=sqrt(sigma2))
  muhat <- mean(x)
  sigmahat <- sd(x)
  SRhat <- (muhat-Rf)/sigmahat
  v[i] <- SRhat
  }
truehist(v)
abline(v=SR,lwd=4) # Ergänzen des wahren Werts



# TEIL III: Simulation der Überdeckungswahrscheinlichkeit der Konfidenzintervalle

v <- matrix(rep(NA,20000),10000,2) # in dieser (10000,2)-Matrix werden die Grenzen der Konfidenzintervalle gespeichert
for(i in 1:10000)
  {
  x <- rnorm(n,mean=mu,sd=sqrt(sigma2))
  muhat <- mean(x)
  sigmahat <- sd(x)
  SRhat <- (muhat-Rf)/sigmahat
  DSR <- matrix(c(1/sigmahat,-0.5*(muhat-Rf)/(sigmahat^3)),1,2)
  Sigma <- matrix(c(sigmahat^2,0,0,2*sigmahat^4),2,2)
  VSR <- DSR%*%Sigma%*%t(DSR)
  v[i,1] <- SRhat-1.96*sqrt(VSR/n) 
  v[i,2] <- SRhat+1.96*sqrt(VSR/n)
  }
print(sum(v[,1]<SR & v[,2]>SR)/10000) # Überdeckungswahrscheinlichkeit (sollte ~0.95 sein)
