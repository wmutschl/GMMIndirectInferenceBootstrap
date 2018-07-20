# Moment estimators for uniform distribution
a <- 0
b <- 1
R <- 10000
ahat <- rep(NA,R)
bhat <- rep(NA,R)
minx <- rep(NA,R)
for(r in 1:R) {
  x <- runif(n=40,a,b)
  ahat[r] <- mean(x)-sqrt(3*var(x))
  bhat[r] <- mean(x)+sqrt(3*var(x))
  minx[r] <- min(x)
  }
par(mfrow=c(2,1))
truehist(ahat)
g <- seq(min(ahat),max(ahat),length=300)
lines(g,dnorm(g,mean(ahat),sd(ahat)))
truehist(bhat)
g <- seq(min(bhat),max(bhat),length=300)
lines(g,dnorm(g,mean(bhat),sd(bhat)))
print("Proportion of impossible estimates:")
print(sum(minx<ahat)/R)
