# Verteilung des MLE, wenn einige Regularitätsannahmen verletzt sind

library(MASS)

for(n in seq(50,1000,by=50)) {
  v <- rep(NA,10000)
  for(i in 1:10000) {
    x <- runif(n,min=0,max=1)
    v[i] <- max(x)
  }
  truehist(v,xlim=c(0.9,1),ylim=c(0,1000),main=paste("n=",n))
}
