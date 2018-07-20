# Nonparametric LM bootstrap test for independence

# Input (or generate) original sample
library(MASS)
n <- 50
z <- mvrnorm(n,c(0,0),matrix(c(1,0.3,0.3,1),2,2))
x <- z[,1]
y <- z[,2]

# Test statistic
Tstat <- cor(x,y)

R <- 10000
Tsharp <- rep(NA,R)

for(r in 1:R) {
  # Resample under H0 (i.e. under independence)
  xx <- sample(x,n,replace=TRUE)
  yy <- sample(y,n,replace=TRUE)
  Tsharp[r] <- cor(xx,yy)
}

Tsharp <- sort(Tsharp)
critlow <- Tsharp[0.025*R]
crithigh <- Tsharp[0.975*R]
if(Tstat<critlow | Tstat>crithigh) print("Reject") else print("Do not reject")
