# Graphical illustration of the three classical tests

# Generate a random draw from an exponential distribution
set.seed(123)
x <- rexp(20,rate=5)
plot(x,rep(0,20),xlab="X",ylab="",main="Sample values (n=20)",axes=F)
axis(1)
box()
v <- locator(1)

# Definition of the log-likelihood function and the score function
loglikelihood <- function(lambda,x) return(length(x)*log(lambda)-lambda*sum(x))
scorefunction <- function(lambda,x) return(length(x)/lambda-sum(x))
rfct <- function(lambda) return(lambda-4)

# Compute functions
lbd <- seq(0.01,10,length=250)
lnL <- rep(NA,length(lbd))
sc <- rep(NA,length(lbd))
rf <- rep(NA,length(lbd))
for(i in 1:length(lbd)) {
  lnL[i] <- loglikelihood(lbd[i],x)
  sc[i] <- scorefunction(lbd[i],x)
  rf[i] <- rfct(lbd[i])
  }

# Wald test
plot(lbd,lnL,type="l",xlab="lambda",ylab="",main="Test idea of Wald test",ylim=c(-5,20))
lines(lbd,rf)
abline(h=0)
text(6.7,17.4,"log-likelihood lnL",pos=4)
text(8.3,3.5,"restriction r",pos=4)
v <- locator(1)
lambdahat <- 1/mean(x)
lines(c(lambdahat,lambdahat),c(0,loglikelihood(lambdahat,x)))
text(lambdahat,-1.5,expression(hat(lambda)[ML]))
v <- locator(1)
lines(c(lambdahat,lambdahat),c(0,rfct(lambdahat)),col="red",lwd=5)
v <- locator(1)

# LR test
plot(lbd,lnL,type="l",xlab="lambda",ylab="",main="Test idea of LR test",ylim=c(-5,20))
lines(lbd,gf)
abline(h=0)
text(6.7,17.4,"log-likelihood lnL",pos=4)
text(8.3,3.5,"restriction r",pos=4)
lines(c(lambdahat,lambdahat),c(0,loglikelihood(lambdahat,x)))
text(lambdahat,-1.5,expression(hat(lambda)[ML]))
text(4,-1.5,expression(hat(lambda)[R]))
v <- locator(1)
lines(c(4,4),c(0,loglikelihood(4,x)))
v <- locator(1)
points(lambdahat,loglikelihood(lambdahat,x),col="red",pch=19)
points(4,loglikelihood(4,x),col="red",pch=19)
v <- locator(1)
lines(c(4,lambdahat),c(loglikelihood(4,x),loglikelihood(4,x)))
lines(c(lambdahat,lambdahat),c(loglikelihood(4,x),loglikelihood(lambdahat,x)),col="red",lwd=5)
v <- locator(1)

# LM test
plot(lbd,lnL,type="l",xlab="lambda",ylab="",main="Test idea LR test",ylim=c(-5,20),col="dark grey")
lines(lbd,rf)
lines(lbd,sc)
abline(h=0)
text(8.3,3.5,"restriction r",pos=4)
text(4,-1.5,expression(hat(lambda)[R]))
text(0.85,19,"score function",pos=4)
v <- locator(1)
lines(c(4,4),c(0,scorefunction(4,x)),col="red",lwd=5)

