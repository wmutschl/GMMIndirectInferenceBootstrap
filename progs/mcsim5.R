# Monte-Carlo-Simulation 5: Verteilung der Dickey-Fuller-Teststatistik
library(MASS)

# Parameter festlegen
n <- 100
R <- 10000
v <- rep(NA,R)

for(r in 1:R)
  {
  # Simulation eines Random-Walk als kumulierte Summe von N(0,1)-ZV
  y <- cumsum(rnorm(n)) 
  deltay <- diff(y)
  lagy <- y[1:(n-1)]
  
  # Einfache Regression ohne Achsenabschnitt
  betahat <- (deltay%*%lagy)/(lagy%*%lagy)
  # Anmerkung: Der Operator %*% berechnet das innere Produkt. Bei Matrizen ist das
  # die Matrix-Multiplikation. Wenn das innere Produkt zweier Vektoren a und b 
  # berechnet wird, ist damit die Summe aller a[i]*b[i] gemeint.

  # Residuen
  uhat <- deltay-betahat*lagy 
  # Geschätzte Störtermvarianz
  sigma2hat <- (uhat%*%uhat)/(n-1) 
  # Standardfehler von betahat, entspricht sigma2hat*(X'X)^(-1)
  SEbetahat <- sqrt(sigma2hat/(lagy%*%lagy)) 
  # t-Teststatistik
  tstat <- betahat/SEbetahat 
  v[r] <- tstat  
  }
  
# Verteilungsfunktion der tau-Statistik
plot(sort(v),seq(0,1,length=R),type="l")
# Ergänzen der Verteilungsfunktion der t-Verteilung
d <- seq(min(v),max(v),length=300)
lines(d,pt(d,df=n-1),col="red")

# Man erkennt: Die t-Verteilung passt für diesen t-Test nicht!
# Wenn in die Regression auch noch ein Achsenabschnitt und/oder
# ein Zeittrend aufgenommen wird, passt die t-Verteilung sogar
# noch deutlich schlechter.
