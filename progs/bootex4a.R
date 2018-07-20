# Parametric Wald bootstrap test for lambda of an exponential distribution

# Inpute (or generate) sample
lambda <- 1
n <- 8
x <- rexp(n,rate=lambda)
# Hypothetical value
lambda0 <- 2
# Test statistic
lambdahat <- 1/mean(x)
Tstat <- lambdahat-lambda0

R <- 10000
Tstar <- rep(NA,R)

for(r in 1:R) {
  # Draw resample
  xx <- rexp(n,rate=lambdahat)
  Tstar[r] <- 1/mean(xx)-lambdahat
}

# Sort Tstar and get confidence intervalls
Tstar <- sort(Tstar)
critlow <- Tstar[0.025*R]
crithigh <- Tstar[0.975*R]
if(Tstat<critlow | Tstat>crithigh)  print("Reject") else print("Do not reject")
