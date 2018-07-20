# Parametrisches 0.95-Perzentil-t-Konfidenzintervall

# Eingabe (oder Generierung) der Originaldaten

lambda <- 1
n <- 8
x <- rexp(n,rate=lambda)
lambdahat <- 1/mean(x)
sigmahat <- lambdahat

# Schleife vorbereiten

B <- 10000
taustar <- rep(NA,B)

for(b in 1:B)
  {
 
  # Schritt 1: Ziehung eines Resamples
  xx <- rexp(n,rate=lambdahat)
  # Schritt 2: lambdastar berechnen und abspeichern
  lambdastar <- 1/mean(xx)
  sigmastar <- lambdastar
  taustar[b] <- sqrt(n)*(lambdastar-lambdahat)/sigmastar
  }

# Sortieren der lambdastar-Werte
taustar <- sort(taustar)

untergrenze <- lambdahat-taustar[0.975*B]*sigmahat/sqrt(n)
obergrenze <- lambdahat-taustar[0.025*B]*sigmahat/sqrt(n)

print(paste("Untergrenze = ",untergrenze))
print(paste("Obergrenze  = ",obergrenze))
