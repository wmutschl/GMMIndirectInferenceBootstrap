# Nonparametric LM bootstrap test for equality of two expectations

# Input (or generate) original sample

n_x <- 15
n_y <- 10
mux <- 2
sdx <- 1
muy <- 2

x <- rnorm(n_x,mean=mux,sd=sdx)
y <- rexp(n_y,rate=1/muy)

# Estimate joint expectation
muhat <- mean(c(x,y))
ax <- mean(x)-muhat
ay <- mean(y)-muhat

R <- 10000
Tsharp <- rep(NA,R)
for(r in 1:R) {
  # Resample under the null hypothesis
  xx <- sample(x-ax,n_x,replace=TRUE)
  yy <- sample(y-ay,n_y,replace=TRUE)
  Tsharp[r] <- (mean(xx)-mean(yy))
}

Tsharp <- sort(Tsharp)
critlow <- Tsharp[0.025*R]
crithigh <- Tsharp[0.975*R]
if(Tstat<critlow | Tstat>crithigh) print("Reject") else print("Do not reject")
