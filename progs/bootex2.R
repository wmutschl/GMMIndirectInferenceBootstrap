# Parametric bootstrap of the bias of lambdahat of an exponential distribution

# Input (or generate) the original sample
lambda <- 2
n <- 800
x <- rexp(n,rate=lambda)
lambdahat <- 1/mean(x)

R <- 1000
lambdastar <- rep(NA,R)

for(r in 1:R) {
  
  # Draw a resample
  xx <- rexp(n,rate=lambdahat)

  # Compute and save lambdastar
  lambdastar[r] <- 1/mean(xx)
}
  
# Calculate the bias
print(mean(lambdastar)-lambdahat)
# Theoertical one can show that Y=sum(X) is Gamma Distributed,
# such that we get a formula for the bias
print(lambda*n/((n-1))-lambda)