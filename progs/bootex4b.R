# Parametric LM bootstrap test for lambda of an exponential distribution
# Inpute (or generate) sample
lambda <- 1
n <- 30
x <- rexp(n,rate=lambda)
# Hypothetical value
lambda0 <- 2
# Test statistic
lambdahat <- 1/mean(x)
Tstat <- lambdahat-lambda0

R <- 10000
Tsharp <- rep(NA,R)

for(r in 1:R) {
  # Draw a resample under H0
  xx <- rexp(n,rate=lambda0)
  Tsharp[r] <- 1/mean(xx)-lambda0
}

# Sort Tsharp and get confidence intervals
Tsharp <- sort(Tsharp)
critlow <- Tsharp[0.025*R]
crithigh <- Tsharp[0.975*R]
if(Tstat<critlow | Tstat>crithigh) print("Reject H0") else print("Do not reject H0")
