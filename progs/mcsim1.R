# Monte-Carlo-Simulation 1: Die Verteilung von lambda-hat einer Exponentialverteilung

# Festlegung einiger Parameter
lambda <- 1
n <- 8
R <- 10000
v <- rep(NA,R) # im Vektor v werden die Schätzwerte abgespeichert

# Programmierung der Schleife
for(r in 1:R)
  {
  # Ziehung der Stichprobe
  x <- rexp(n,rate=lambda)
  # Berechnung und Speicherung des Schätzwerts
  v[r] <- 1/mean(x)
  }

# Erwartungswert
print(mean(v))

# Varianz
print(var(v))

# Dichte
plot(density(v))

# Verteilungsfunktion
plot(sort(v),seq(0,1,length=length(v)),type="l")
lines(sort(v),pnorm(sort(v),mean=1,sd=sqrt(1/(n*lambda^2))),col="red")

# 0.95-Intervall
print(quantile(v,c(0.025,0.975)))
