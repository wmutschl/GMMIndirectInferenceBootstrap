# Empirische Illustration:
# Hat sich die Volatilität des DAX geändert?
library(MASS)
library(tseries)

# DAX-Kurse einlesen (Spalte1: Datum, Spalte2: Index)
setwd("c:/skripte/inferenz/progs")
x <- read.csv2("dax.csv")

# Renditen berechnen und an die Datenmatrix anhängen
r <- c(NA,diff(log(x[,2])))
x <- cbind(x,r)

# Feiertage löschen (d.h wenn Rendite=0)
x <- x[r!=0,]

# Nur die letzten 500 Tage betrachten
x <- x[(dim(x)[1]-499):dim(x)[1],]
r <- x[,3]

#####################################
# I. Schätzung des GARCH(1,1)-Modells
#####################################

muhat <- mean(r)
mod <- garch(r-muhat)
alpha0hat <- coefficients(mod)[1]
alpha1hat <- coefficients(mod)[2]
beta1hat <- coefficients(mod)[3]

# Schätzung der bedingten Varianzen
sigma2hat <- rep(NA,500)
sigma2hat[1] <- alpha0hat/(1-alpha1hat-beta1hat)
for(i in 2:500)
  sigma2hat[i] <- alpha0hat+alpha1hat*(r[i-1]-muhat)^2+beta1hat*sigma2hat[i-1]
# Alternativ: sigma2hat <- mod$fitted.values[,1]^2

# Residuen berechnen
epshat <- (r-muhat)/sqrt(sigma2hat)
# Alternativ: epshat <- residuals(mod)

# Standardisierung der Residuen
epshat <- (epshat-mean(epshat))/sd(epshat)

####################
# II. Bootstrap-Test
####################

# Berechnung der Original-Teststatistik
F <- var(r[1:250])/var(r[251:500])

# Bootstrap-Schleife
B <- 1000
Fstar <- rep(NA,B)
for(b in 1:B)
  {
  # Resample aus den Residuen
  epsstar <- c(NA,sample(epshat,499,replace=TRUE))

  # Generiere Resample-Rendite-Verlauf
  xstar <- rep(NA,500)
  sigma2star <- rep(NA,500)
  rstar <- rep(NA,500)

  xstar[1] <- r[1]-muhat
  sigma2star[1] <- sigma2hat[1]
  rstar[1] <- r[1]

  for(i in 2:500)
    {
    sigma2star[i] <- alpha0hat+alpha1hat*xstar[i-1]^2+beta1hat*sigma2star[i-1]
    xstar[i] <- sqrt(sigma2star[i])*epsstar[i]
    rstar[i] <- xstar[i]+muhat
    }
  
  # Bootstrap-Teststatistik
  Fstar[b] <- var(rstar[1:250])/var(rstar[251:500])
  }

Fstar <- sort(Fstar)

print(paste("Wert der Teststatistik  = ",F))
print(paste("Unterer kritischer Wert = ",Fstar[0.025*B]))
print(paste("Oberer kritischer Wert  = ",Fstar[0.975*B]))
