# Bootstrap-Konfidenzintervall für den Autokorrelationskoeffizienten
# erster Ordnung eines stationären stochastischen Prozesses

# Daten einlesen (bzw. generieren)
mu <- 5
a1 <- 0.7
sigma2 <- 1
T <- 241
x <- filter(rnorm(T,0,sd=sqrt(sigma2)),a1,method="recursive",init=rnorm(1,0,sd=sqrt(sigma2/(1-a1^2))))+mu
y <- cbind(x[1:(T-1)],x[2:T])

# Schätzung
muy1hat <- mean(y[,1])
a1hat <- sum((y[,1]-muy1hat)*(y[,2]-muy1hat))/sum((y[,1]-muy1hat)^2)

# MBB vorbereiten
L <- 12
B <- 1000
a1star <- rep(NA,B)
ystar <- y

for(b in 1:B)
  {
  indices <- sample(1:(T-L),(T-1)/L,replace=TRUE)
  for(v in 0:(L-1))
    ystar[seq(1+v,T-L+1+v,by=L),] <- y[indices+v,]
  muy1star <- mean(ystar[,1])
  a1star[b] <- sum((ystar[,1]-muy1star)*(ystar[,2]-muy1star))/sum((ystar[,1]-muy1star)^2)
  }

a1star <- sort(a1star)

print(paste("Punktschätz = ",a1hat))
print(paste("Untergrenze = ",2*a1hat-a1star[0.975*B]))
print(paste("Obergrenze  = ",2*a1hat-a1star[0.025*B]))
