# Monte-Carlo-Simulation 2: Die Verteilung des IV-Schätzers betahhatIV in kleinen Stichproben
library(MASS)

# Festlegung einiger Parameter
n <- 25
R <- 1000
beta0 <- matrix(c(10,3))
Mean <- c(0,10,1,2,3,4,5)
Sigma <- matrix(c(1,0.6,0,0,0,0,0,0.6,1,0.7,0.6,0.5,0.4,0.3,0,0.7,1,0.7,0.6,0.5,0.4,
                  0,0.6,0.7,1,0.7,0.6,0.5,0,0.5,0.6,0.7,1,0.7,0.6,0,0.4,0.5,0.6,0.7,
                  1,0.7,0,0.3,0.4,0.5,0.6,0.7,1),7,7)


############################################################################
# I. Schätzung des Erwartungswert und der Verteilungsfunktion von betahatOLS
############################################################################

# In der Matrix v werden die Schätzwerte abgespeichert (Achsenabschnitt und Steigung)
v <- matrix(NA,R,2)
for(r in 1:R)
  {
  # Ziehung aus der multivariaten Normalverteilung und Zuweisungen zu Variablen
  # (hier wird eine 7x7-Matrix gezogen, obwohl nur die ersten zwei Spalten 
  # gebraucht werden, das könnte man verschlanken)
  M <- mvrnorm(n,Mean,Sigma)
  u <- M[,1,drop=FALSE]           # Durch die drop-Option wird erreicht, dass class(u)="matrix" ist
  X <- cbind(1,M[,2,drop=FALSE])  # (normalerweise wird eine einspaltige Matrix zu einem Vektor konvertiert)
  
  # Berechnung der endogenen Variablen y
  y <- X%*%beta0+u
  
  # Berechnung des OLS-Schätzers
  betahatOLS <- solve(t(X)%*%X)%*%t(X)%*%y
  v[r,] <- betahatOLS
  }
# Erwartungswert von betahatOLS_1
print(mean(v[,1]))
# Erwartungswert von betahatOLS_2
print(mean(v[,2]))
# Beide Schätzer sind offenbar verzerrt
# (die wahren Werte sind 10 und 3)


###############################################
# II. Schätzung für ein einzelnes Instrument Z1
###############################################

# Zunächst exakt wie unter I.
v <- matrix(NA,R,2)
for(r in 1:R)
  {
  M <- mvrnorm(n,Mean,Sigma)
  u <- M[,1,drop=FALSE]         
  X <- cbind(1,M[,2,drop=FALSE])
  Z1 <- cbind(1,M[,3,drop=FALSE]) # Zu den Instrumentvariablen gehört auch der Achsenabschnitt!
  y <- X%*%beta0+u
  P <- Z1%*%solve(t(Z1)%*%Z1)%*%t(Z1)
  betahatIV1 <- solve(t(X)%*%P%*%X)%*%t(X)%*%P%*%y
  v[r,] <- betahatIV1
  }
# Erwartungswert von betahatIV_1
print(mean(v[,1]))
# Erwartungswert von betahatIV_2
print(mean(v[,2]))
# Verteilungsfunktion von betahatIV_2
plot(sort(v[,2]),seq(0,1,length=R),type="l")
abline(v=beta0[2],lwd=3) # Ergänzen des wahren Werts ("v=" steht für vertikale Linie)


#############################################
# III. Schätzung für alle fünf Instrumente Z5
#############################################

# Zunächst exakt wie unter II.
v <- matrix(NA,R,2)
for(r in 1:R)
  {
  M <- mvrnorm(n,Mean,Sigma)
  u <- M[,1,drop=FALSE]         
  X <- cbind(1,M[,2,drop=FALSE])
  Z5 <- cbind(1,M[,3:7])
  y <- X%*%beta0+u
  P <- Z5%*%solve(t(Z5)%*%Z5)%*%t(Z5)
  betahatIV5 <- solve(t(X)%*%P%*%X)%*%t(X)%*%P%*%y
  v[r,] <- betahatIV5
  }
# Erwartungswert von betahatIV_1
print(mean(v[,1]))
# Erwartungswert von betahatIV_2
print(mean(v[,2]))
# Verteilungsfunktion von betahatIV_2 in die Grafik von II. ergänzen
lines(sort(v[,2]),seq(0,1,length=R),type="l",col="red")
# Man erkennt: Die fünf Instrumente bringen nicht unbedingt eine Verbesserung.
# Der Bias steigt wieder an. Zur Erläuterung s. z.B. Davidson, MacKinnon (1993),
# Estimation and Inference in Econometrics.
