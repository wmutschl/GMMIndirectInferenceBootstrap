# Monte-Carlo-Simulation 4: GMM-Schätzer in kleinen Stichproben

# Parameter festlegen
N <- 10
T <- 50
R <- 1000

v <- matrix(NA,R,2)
for(r in 1:R)
  {
  y <- matrix((rlnorm(N*T,0,1)-1.648721)/2.161197,N,T)
  sigma2i <- apply(y,1,var)
  X <- matrix(1,N,1) 

  # Schätzung ohne Gewichtung (hier könnte man auch einfach schreiben: thetahat1 <- mean(sigma2i) 
  thetahat1 <- solve(t(X)%*%X)%*%t(X)%*%sigma2i
  v[r,1] <- thetahat1

  # Schätzung mit asymptotisch optimaler GMM-Gewichtung
  # Zuerst Schätzung der Kovarianzmatrix Omegahat
  Omegahat <- matrix(NA,N,N)
  for(i in 1:N)
    {
    for(j in i:N)
      {
      Omegahat[i,j] <- sum((y[i,]-mean(y[i,]))^2*(y[j,]-mean(y[j,]))^2)/(T-1)-sigma2i[i]*sigma2i[j]
      Omegahat[j,i] <- Omegahat[i,j] # wegen Symmetrie
      }
    }
  # Berechnung der Inversen
  invOmega <- solve(Omegahat)
  thetahat2 <- solve(t(X)%*%invOmega%*%X)%*%t(X)%*%invOmega%*%sigma2i
  v[r,2] <- thetahat2
  }

# Verteilungsfunktion von thetahat1
plot(sort(v[,1]),seq(0,1,length=R),type="l")
# Einfügen der Verteilungsfunktion von thetahat2
lines(sort(v[,2]),seq(0,1,length=R),col="red")
# Einfügen des wahren Werts
abline(v=1,lwd=3)

# Berechnen der Mean-Squared-Errors
print(mean((v[,1]-1)^2))
print(mean((v[,2]-1)^2))

# Man erkennt, dass der GMM-Schätzer in dieser Situation total versagt.
# Die Ergebnisse entsprechen recht gut den in Tabelle 1 (Seite 40)
# aufgeführten Ergebnissen von Altonji und Segal (Working-Paper-Version).
# Sie finden das Working Paper als pdf auf der Internetseite der Vorlesung.
