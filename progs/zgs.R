# R-Programm zur Illustration des zentralen Grenzwertsatzes

# Aufrufen des Pakets "MASS", damit einige nützliche weitere
# Befehle zur Verfügung stehen (insb. truehist für Histogramme)

library(MASS)

# Die folgenden beiden Parameter dienen der Standardierung
# des Durchschnitts, sie hängen mit dem rNAME-Befehl weiter
# unten zusammen.

mu <- 1
sigma <- 1

R <- 10000 # Anzahl der zu ziehenden Stichproben (jeweils der Länge n)
g <- seq(-5,5,length=500) # Gitterdefinition

# Die Schleife läuft über die Stichprobenumfänge
# (die R-fache Ziehung der Stichproben wird matriziell programmiert)

for(n in 1:50)
  {

  # Die Matrix x hat R Zeilen und n Spalten, jede Zeile ist
  # eine eigene (unabhängige) Stichprobe; die Verteilung sollte
  # nicht die Normalverteilung sein (sonst ist es langweilig)

  x <- matrix(rexp(n*R),R,n)

  # Der apply-Befehl berechnet für jede Zeile von x den
  # Durchschnitt und gibt alle R Durchschnitte als Vektor zurück

  xbar <- apply(x,1,mean)
  z <- sqrt(n)*(xbar-mu)/sigma # Standardisierung
  truehist(z) # Erzeugen einer Grafik mit dem Histogramm der z-Werte
  lines(g,dnorm(g)) # Einfügen der N(0,1)-Dichte
  }
  
# Mögliche Verbesserungen:
# (1) Erweitern des truehist-Befehls:
#     truehist(z,xlim=c(-5,5), ylim=c(0,0.5), main=paste("n = ",n))
# (2) Erweitern des lines-Befehls:
#     lines(g,dnorm(g),col="red",lwd=3)

# Andere Verteilungen:
# (1) Recheckverteilung:
#     mu <- 0.5
#     sigma <- sqrt(1/12)
#     x <- matrix(runif(n*R),R,n)
#
# (2) Student-t-Verteilung mit 3 df: 
#     mu <- 0
#     sigma <- sqrt(3)
#     x <- matrix(rt(n*R,df=3),R,n)
