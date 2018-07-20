# Parametric bootstrap 0.95-confidence interval 
# for lambdahat of an exponential distribution

# Input (or generate) the original sample
lambda <- 1
n <- 8
x <- rexp(n,rate=lambda)
lambdahat <- 1/mean(x)

R <- 10000
lambdastar <- rep(NA,R)
for(r in 1:R) {
  # Draw a resample
  xx <- rexp(n,rate=lambdahat)
  # Compute and save lambdastar
  lambdastar[r] <- 1/mean(xx)
}
# Sort lambdastar
lambdastar <- sort(lambdastar)
Lower_Boot <- 2*lambdahat-lambdastar[0.975*R]
Upper_Boot <- 2*lambdahat-lambdastar[0.025*R]

#Compare with normal approximated and exact confidence intervalls
Lower_Approx <- lambdahat*(1-1.96/sqrt(n))
Upper_Approx <- lambdahat*(1+1.96/sqrt(n))
Lower_Exact <- lambdahat*qchisq(0.025,2*n)/(2*n)
Upper_Exact <- lambdahat*qchisq(0.975,2*n)/(2*n)
CI <- matrix(c(Lower_Approx,Lower_Boot,Lower_Exact,Upper_Approx,Upper_Boot,Upper_Exact),3,2)
rownames(CI) <- c("Approx","Boot","Exact")
colnames(CI) <- c("Lower","Upper")
print(CI)