# Simultanes Konfidenzband

# Wahre Regressionsgerade

a <- 0
b <- 1

n <- 20
x <- seq(0,1,length=n)
u <- rnorm(n)
y <- a+b*x+u
plot(x,y,ylim=c(-2,2.5))
abline(a,b,col="light grey")

# Definition der Gitterpunkte
x0 <- seq(0,1,length=100)

# Parameter schätzen
mod <- lm(y~x)
ahat <- coefficients(mod)[1]
bhat <- coefficients(mod)[2]
uhat <- residuals(mod)
sigma2hat <- sum(uhat^2)/(n-2)
f <- sqrt(1/n+(x0-mean(x))^2/sum((x-mean(x))^2))
abline(ahat,bhat)

B <- 1000
Deltastar <- matrix(NA,B,length(x0))

# Erzeuge die Bootstrap-Kurven
for(j in 1:B)
  {
  # Resample aus den Residuen
  ustar <- sample(uhat,n,replace=TRUE)
  ystar <- ahat+bhat*x+ustar
  # Bootstrap: Parameter schätzen
  bmod <- lm(ystar~x)
  astar <- coefficients(bmod)[1]
  bstar <- coefficients(bmod)[2]
  ustar <- residuals(bmod)
  sigma2star <- sum(ustar^2)/(n-2)

  Deltastar[j,] <- ((astar+bstar*x0)-(ahat+bhat*x0))/sqrt(sigma2star)
  }
  
# Suche die Schätzer c1star und c2star
c1star <- 1.96
while(sum(apply(Deltastar,1,function(z) all(z < c1star*f)))/B<0.975) c1star<-c1star+0.005
c2star <- 1.96
while(sum(apply(Deltastar,1,function(z) all(z > -c2star*f)))/B<0.975) c2star<-c2star+0.005

# Obere und untere Konfidenzgrenze
confbandlow <- ahat+bhat*x0-c1star*sqrt(sigma2hat)*f
confbandhigh <- ahat+bhat*x0+c2star*sqrt(sigma2hat)*f

# Zeichne die Regressionsgerade und das Konfidenzband
lines(x0,confbandhigh)
lines(x0,confbandlow)
