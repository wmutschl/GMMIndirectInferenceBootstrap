# Why are naive bootstrap confidence intervals often wrong?

windows(12,8)
par(yaxs="i",xpd=TRUE,cex=1.2)
d <- seq(0,30,length=500)
mue <- 1.6
sigm <- 0.5
shift <- 5
muelog <- exp(1.6+0.5^2/2)
estimated <- muelog +shift

plot(curve(dlnorm(x,mue,sigm),from=0,to=30),type="l",xlab="estimator",ylab="Density of estimator",ylim=c(0,0.2),main="Naive bootstrap confidence intervals")
text(1.3,0.07,"density",pos=2)
lines(rep(muelog,2),c(0,0.19),lwd=2)
text(5.5,0.17,"true value",pos=4)
locator(1)
lines(rep(estimated,2),c(0,0.19),col="green",lwd=2)
text(5.5+shift,0.17,"estimated value",col="green",pos=4)
locator(1)
curve(dlnorm(x-shift,mue,sigm),from=0,to=30,add=T,lty="dashed",col="green",lwd=2)
#lines(d,dlnorm(d-5,1.6,0.5),lty="dashed",col="green",lwd=2)
text(1.3+2.25*shift,0.07,"bootstrap density",col="green",pos=4)


### Stuff
locator(1)
plot(curve(dlnorm(x,mue,sigm),from=0,to=30),type="l",xlab="estimator",ylab="Density of estimator",ylim=c(0,0.2),main="Percentile bootstrap confidence intervals")
lines(rep(muelog,2),c(-0.005,0.005),lwd=5)
lines(rep(estimated,2),c(-0.005,0.005),lwd=5,col="green")
curve(dlnorm(x-shift,mue,sigm),from=0,to=30,add=T,lty="dashed",col="green",lwd=2)

locator(1)
curve(dlnorm(2*estimated-x-shift,mue,sigm),from=0,to=30,add=T,lty="solid",col="dark green",lwd=2)
lines(rep(estimated,2),c(0,0.2),col="green")
