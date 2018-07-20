# Punktweises Konfidenzband

# Wahre Regressionsgerade

a <- 0
b <- 1

n <- 20
x <- seq(0,1,length=n)


# I. Punktweises Prognoseband zeichnen
x0 <- seq(0,1,length=100)
u <- rnorm(n)
y <- a+b*x+u
mod <- lm(y~x)
ahat <- coefficients(mod)[1]
bhat <- coefficients(mod)[2]
uhat <- residuals(mod)
sigma2hat <- sum(uhat^2)/(n-2)
se <- sqrt(sigma2hat*(1/n+(x0-mean(x))^2/sum((x-mean(x))^2))) 
thetalow <- ahat+bhat*x0-1.96*se
thetahigh <- ahat+bhat*x0+1.96*se
plot(x0,thetalow,type="l",ylim=range(c(thetalow,thetahigh)))
lines(x0,thetahigh)
abline(ahat,bhat,col="red")

abline(a,b,col="light grey")
points(x,y)


# II. EINE Prognosestelle festlegen und die Überdeckungswahrscheinlichkeit auswerten
x0 <- 0.06

N <- 1000
S <- 0
for(i in 1:N)
  {

  # y-Werte generieren
  u <- rnorm(n)
  y <- a+b*x+u

  # Regressionsgerade schätzen
  mod <- lm(y~x)
  ahat <- coefficients(mod)[1]
  bhat <- coefficients(mod)[2]
  uhat <- residuals(mod)
  sigma2hat <- sum(uhat^2)/(n-2)
  se <- sqrt(sigma2hat*(1/n+(x0-mean(x))^2/sum((x-mean(x))^2))) 

  thetalow <- ahat+bhat*x0-1.96*se
  thetahigh <- ahat+bhat*x0+1.96*se

  if(thetalow<(a+b*x0) & thetahigh>(a+b*x0)) S <- S+1
  }
  
print(paste("Die Überdeckungswahrscheinlichkeit an einer Stelle x0 ist ",S/N))


# III. MEHRERE Prognosestellen festlegen und die Überdeckungswahrscheinlichkeit auswerten
x0 <- seq(0,1,length=100)

S <- 0
for(i in 1:N)
  {

  # y-Werte generieren
  u <- rnorm(n)
  y <- a+b*x+u

  # Regressionsgerade schätzen
  mod <- lm(y~x)
  ahat <- coefficients(mod)[1]
  bhat <- coefficients(mod)[2]
  uhat <- residuals(mod)
  sigma2hat <- sum(uhat^2)/(n-2)
  se <- sqrt(sigma2hat*(1/n+(x0-mean(x))^2/sum((x-mean(x))^2))) 

  thetalow <- ahat+bhat*x0-1.96*se
  thetahigh <- ahat+bhat*x0+1.96*se

  if(all(thetalow<(a+b*x0)) & all(thetahigh>(a+b*x0))) S <- S+1
  }
  
print(paste("Die Überdeckungswahrscheinlichkeit an mehreren Stellen x0 ist ",S/N))
