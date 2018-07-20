# A simple example of indirect inference estimation
# Estimate a MA(1) model with auxiliary model AR(3)
# Comparison of the distribution of betahat_ML and betahat_II


# Function to compute indirect inference estimator of MA(1) with AR(3) as auxiliary model
MA1_indirect <- function(truedata,H,W,start_val,seed_nr){
  # Estimate parameters theta from auxiliary model AR(3) for true data
  theta_hat <- coefficients(arima(truedata,order=c(3,0,0),include.mean=F,method="ML"))
  
  dist_fct <- function(betta) {
    set.seed(seed_nr)
    theta_hatsim <- matrix(NA,3,H)
    for(h in 1:H) {
      simdata <- filter(rnorm(TT+1),filter=c(1,betta),method="c")[-(TT+1)] # simulate new data from MA(1) model
      theta_hatsim[,h] <- coefficients(arima(simdata,order=c(3,0,0),include.mean=F),method="ML") #estimate auxiliary model
    }
    theta_tilde <- rowMeans(theta_hatsim)
    Q <- t(theta_hat - theta_tilde) %*% W %*% (theta_hat - theta_tilde)
    return(Q)
  }
  
  # Compute indirect inference estimate
  beta_II <- optimize(dist_fct,interval=start_val)$minimum
  return(beta_II)
}

# True data
TT <- 250
beta_true <- 0.5
truedata <- filter(rnorm(TT+1),filter=c(1,beta_true),method="c")[-(TT+1)]

# Compute indirect inferece estimate of MA(1) with AR(3) as auxiliary model
MA1_indirect(truedata=truedata,H=100,W=diag(3),start_val = c(0,0.8),seed_nr = as.integer(runif(1)*2e7))
# For comparison: Compute ML estimate of MA(1)
coefficients(arima(truedata,order=c(0,0,1),include.mean=FALSE),method="ML")

# Get distribution of indirect inference estimator
R <- 100
Z <- matrix(NA,R,2)
fb <- txtProgressBar(min=0, max=R, style=3)
for(i in 1:R) {
  # True data sample
  truedata <- filter(rnorm(TT+1),filter=c(1,beta_true),method="c")[-(TT+1)]
  # Compute indirect inferece estimate of MA(1) with AR(3) as auxiliary model
  Z[i,1] <- MA1_indirect(truedata=truedata,H=10,W=diag(3),start_val = c(0,0.8),seed_nr = as.integer(runif(1)*2e7))
  # For comparison: Compute ML estimate of MA(1)
  Z[i,2] <- coefficients(arima(truedata,order=c(0,0,1),include.mean=FALSE),method="ML")
  
  setTxtProgressBar(fb,i)
}
close(fb)


# Compare distributions
print(apply(Z,2,mean))
print(apply(Z,2,sd))
plot(ecdf(Z[,1]),main="IndInf and ML Estimators")
lines(ecdf(Z[,2]),col="red")
legend(locator(1),c("IndInf","ML"),fill=c("black","red"))
library(MASS)
par(mfrow=c(1,2))
truehist(Z[,1],main="IndInf")
truehist(Z[,2],main="ML")

