# Nonparametric Wald bootstrap test for equality of two expectations

# Input (or generate) original sample
n_x <- 15
n_y <- 10
mux <- 2
sdx <- 1
muy <- 2

x <- rnorm(n_x,mean=mux,sd=sdx)
y <- rexp(n_y,rate=1/muy)

# Compute test statistic
Tstat <- (mean(x)-mean(y))/sqrt(var(x)+var(y))

R <- 10000
Tstar <- rep(NA,R)
for(r in 1:R)  {
  # Resample under the alternative hypothesis
  xx <- sample(x,n_x,replace=TRUE)
  yy <- sample(y,n_y,replace=TRUE)
  Tstar[r] <- (mean(xx)-mean(yy))/sqrt(var(xx)+var(yy))
}

Tstar <- sort(Tstar)
critlow <- Tstar[0.025*R]
crithigh <- Tstar[0.975*R]
if(Tstat<critlow | Tstat>crithigh) print("Reject") else print("Do not reject")
