# Nonparametric bootstrap of the standard error of xbar

# Input (or generate) the original sample

n <- 20
mue <- 50
sigm <- 3
x <- rnorm(n,mean=mue,sd=sigm)

R <- 1000
xxbar <- rep(NA,R)

for(r in 1:R) {
  # Draw a resample
  xx <- sample(x,n,replace=TRUE)
  # Compute and save the mean
  xxbar[r] <- mean(xx)
}
  
# Compute the standard error
sqrt(1/(R-1)*sum((xxbar-mean(xxbar))^2))
print(sd(xxbar))

# Comparison with theoretical standard error due to CLT
print(sd(x)/sqrt(n))
print(sigm/sqrt(n))
