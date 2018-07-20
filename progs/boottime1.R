# Bootstrap-Konfidenzintervall für den Erwartungswert eines AR(1)

# Eingabe oder Generieren der Daten
mu <- 5
alpha1 <- 0.9
sigma2eps <- 1
T <- 100

x <- rep(NA,100)
x[1] <- rnorm(1,sd=sqrt(sigma2eps/(1-alpha1^2)))
for(i in 2:T)
  x[i] <- mu+alpha1*(x[i-1]-mu)+rnorm(1,sd=sqrt(sigma2eps))

# Schätzer
muhat <- mean(x)
alpha1hat <- sum((x[-1]-muhat)*(x[-T]-muhat))/sum((x-muhat)^2)
epshat <- (x[-1]-muhat)-alpha1hat*(x[-T]-muhat)

estimate <- arima(x,order=c(1,0,0),include.mean=TRUE)

# Zeitreihen-Bootstrap
B <- 1000
mustar <- rep(NA,B)
for(b in 1:B)
  {
  # Resample ziehen
  epsstar <- c(NA,sample(epshat,T-1,replace=TRUE))
  xstar <- rep(NA,T)
  xstar[1] <- x[1]
  for(i in 2:T)
    xstar[i] <- muhat+alpha1hat*(xstar[i-1]-muhat)+epsstar[i]

  # Bootstrap-Schätzer
  mustar[b] <- mean(xstar)
  }
mustar <- sort(mustar)

print("Zeitreihen-Bootstrap")
print(paste("Wahrer Wert = ",mu))
print(paste("muhat       = ",muhat))
print(paste("Untergrenze = ",2*muhat-mustar[0.975*B]))
print(paste("Obergrenze  = ",2*muhat-mustar[0.025*B]))


# Zum Vergleich: iid-Bootstrap
print("Zum Vergleich: iid-Bootstrap")
B <- 1000
mustar <- rep(NA,B)
for(b in 1:B)
  {
  xstar <- sample(x,T,replace=TRUE)
  mustar[b] <- mean(xstar)
  }
mustar <- sort(mustar)
print(paste("Untergrenze = ",2*muhat-mustar[0.975*B]))
print(paste("Obergrenze  = ",2*muhat-mustar[0.025*B]))
