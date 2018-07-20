# Bootstrap der Residuen

# Eingabe (oder Generieren) der Beobachtungen (x,y)
x <- runif(25)*5
u <- rexp(25)-1
y <- 5+2*x+u
n <- length(x)

# OLS-Schätzung
mod <- lm(y~x)
uhat <- residuals(mod)
ahat <- coefficients(mod)[1]
bhat <- coefficients(mod)[2]
SEbhat <- sqrt(vcov(mod)[2,2])

# Perzentil-t-Bootstrap
R <- 1000
taustar <- rep(NA,R)
for(r in 1:R)
  {
  ustar <- sample(uhat,n,replace=TRUE)
  ystar <- ahat+bhat*x+ustar
  bootmod <- lm(ystar~x)
  bstar <- coefficients(bootmod)[2]
  SEbstar <- sqrt(vcov(bootmod)[2,2])
  taustar[r] <- (bstar-bhat)/SEbstar
}
  
taustar <- sort(taustar)
Lower_Boot <- bhat-taustar[0.975*R]*SEbhat
Upper_Boot <- bhat-taustar[0.025*R]*SEbhat
Lower_Approx <- bhat-1.96*SEbhat
Upper_Approx <- bhat+1.96*SEbhat
CI <- matrix(c(Lower_Approx,Lower_Boot,Upper_Approx,Upper_Boot),2,2)
rownames(CI) <- c("Approx","Boot")
colnames(CI) <- c("Lower","Upper")
print(CI)