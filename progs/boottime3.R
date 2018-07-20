# Bootstrap-Konfidenzintervall für den Erwartungswert eines 
# stationären stochastischen Prozesses


S <- 0
for(i in 1:100)
{


# Daten einlesen (bzw. generieren)
mu <- 5
a1 <- 0.7
sigma2 <- 1
T <- 240
x <- filter(rnorm(T,0,sd=sqrt(sigma2)),a1,method="recursive",init=rnorm(1,0,sd=sqrt(sigma2/(1-a1^2))))+mu

# Schätzung
muhat <- mean(x)

# MBB vorbereiten
L <- 20
B <- 1000
mustar <- rep(NA,B)
xstar <- x

for(b in 1:B)
  {
  indices <- sample(1:(T-L+1),T/L,replace=TRUE)
  for(v in 0:(L-1))
    xstar[seq(1+v,T-L+1+v,by=L)] <- x[indices+v]
  mustar[b] <- mean(xstar)
  }

mustar <- sort(mustar)

#print(paste("Untergrenze = ",2*muhat-mustar[0.975*B]))
#print(paste("Obergrenze  = ",2*muhat-mustar[0.025*B]))

if(2*muhat-mustar[0.975*B]<mu & 2*muhat-mustar[0.025*b]>mu) S <- S+1
}
print(S/100)
