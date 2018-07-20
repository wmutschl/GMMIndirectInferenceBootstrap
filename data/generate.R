# Generate datasets for the course "Statistical Inference"

library(MASS)
setwd("c:/temp")

# expgrowth.csv
alpha0 <- 2
beta0 <- 0.05
x <- runif(100,min=0,max=40)
u <- rnorm(100,sd=3)
y <- exp(alpha0+beta0*x)+u
z <- data.frame(y=round(y,2),x=round(x,2))
write.csv(z,file="expgrowth.csv",row.names=F)

# tobitbsp.csv
sigma <- 1/8.022
beta0 <- c(1.3392,-0.2246,0.0350)*sigma
MM <- matrix(c(1, 3.8095238,0.36763265,0.08360408,
 3.80952381,16.4843537,1.77356463,0.28244626,
 0.36763265,1.7735646,0.89524803,0.02797143,
 0.08360408,0.2824463,0.02797143,0.01784093),4,4)
Sigma <- MM[2:3,2:3]-outer(MM[1,2:3],MM[1,2:3])
XX <- mvrnorm(n=735,MM[1,2:3],Sigma)
y <- beta0[1]+XX%*%beta0[2:3]+rnorm(n=735,sd=sigma)
x2 <- round(XX[,1],0)
x2[x2>6] <- 6
x2[x2==0] <- 1
z <- data.frame(y=round(y,5),x1=rep(1,735),x2=x2,x3=round(XX[,2],5))
z$y[z$y<0] <- 0
write.csv(z,file="tobitbsp.csv",row.names=F)

# censoredln.csv
x <- round(rlnorm(200,2,0.5),4)
x[x>12] <- 12
write.csv(x,file="censoredln.csv",row.names=F)

# mroz.csv
x <- read.csv("mroz.raw",header=F,sep="",na.strings=".")
names(x) <- c("inlf","hours","kidslt6","kidsge6","age","educ","wage","repwage","hushrs","husage","huseduc","huswage","faminc","mtr","motheduc","fatheduc","unem","city","exper","nwifeinc","lwage","expersq")
write.csv(x,file="mroz.csv",row.names=F)

# spells.csv
set.seed(1234)
n <- 1000
X <- cbind(1,runif(n)<0.4,rnorm(n,1,0.5))
alpha <- 1.3
beta0 <- c(1,0.7,0.3)
theta <- exp(X%*%beta0)
tt <- rweibull(n,shape=alpha,scale=1/theta)
tt[tt>0.5] <- 0.5
x <- data.frame(duration=round(tt,5),const=X[,1],X1=X[,2],X2=round(X[,3],5))
write.csv(x,file="spells.csv",row.names=F)

# panelre.csv
set.seed(1234)
NN <- 75
TT <- 10
mu <- 3
beta0 <- c(2,-2)
sigma2a <- 2
sigma2e <- 1
pid <- rep(1:NN,each=TT)
tt <- rep(1:TT,times=NN)
x <- data.frame(pid=pid,tt=tt,y=rep(NA,NN*TT),x1=rep(NA,NN*TT),x2=rep(NA,NN*TT))
for(i in 1:NN) {
  x1 <- mvrnorm(1,runif(10,min=10,max=20),15*(diag(10)+matrix(0.4,10,10)))
  x2 <- as.numeric(runif(10)<0.5)
  x[((i-1)*TT+1):(i*TT),4] <- x1
  x[((i-1)*TT+1):(i*TT),5] <- x2
  alphai <- rnorm(1,0,sqrt(sigma2a))
  eps <- rnorm(TT,0,sqrt(sigma2e))
  y <- mu+beta0[1]*x1+beta0[2]*x2+alphai+eps
  x[((i-1)*TT+1):(i*TT),3] <- y
  }
x <- round(x,2)
write.csv(x,file="panelre.csv",row.names=F)

# gaussian.csv
library(copula)
x <- rcopula(tCopula(0,4,df=1.2),n=10000)
V <- qnorm(x)
x <- data.frame(V1=V[,1],V2=V[,2],V3=V[,3],V4=V[,4])
write.csv(x,file="c:/temp/gaussian.csv",row.names=F)

# DMacK1.csv
library(copula)
set.seed(1234)
x <- rcopula(frankCopula(4,2),n=250)
x[,1] <- qnorm(x[,1])
x[,2] <- qexp(x[,2])
b1 <- 1
b2 <- 3
y <- b1+b2*x[,1]+x[,2]/b2+rnorm(250,sd=1)
x1 <- x[,1]
x2 <- x[,2]
write.csv(round(cbind(y,x1,x2),5),file="DMacK1.csv",row.names=F)

# AR(3)
set.seed(1234)
b <- c(0.8,-0.5,0.2)
x <- filter(rnorm(403),filter=b,method="r",init=c(rnorm(3)))
x <- data.frame(time=-2:400,x=round(x,5))
write.csv(x,file="ar3bsp.csv",row.names=F)

# censoredln.csv
set.seed(1234)
y <- rlnorm(n=500,2.1,0.5)
y[y>=12] <- 12
write.csv(round(y,5),file="censoredln.csv",row.names=F)

# players.csv
set.seed(1234)
n <- 1000
position <- runif(n)<0.3
age <- round(runif(n,min=10,max=15)+runif(n,min=7,max=16),0)
age2 <- age^2
training <- round(rnorm(n,mean=15,sd=2),0)
salary <- round(runif(n,min=100,max=1000),1)
bonus <- round(ifelse(position==0,rpois(n,1)*10,round(runif(n,min=5,max=50))),1)
indx <- -2.4+0.5*position+0.095*age-0.00167*age2+0.02*training+0.001*salary+0.02*bonus
mu <- exp(indx)
y <- rep(NA,n)
for(i in 1:n) y[i] <- rpois(1,mu[i])
X <- data.frame(goals=y,position=as.numeric(position),age,training,salary,bonus)
write.csv(X,file="players.csv",row.names=F)

# arch1bsp.csv
set.seed(1234)
omega <- 0.2
alpha <- 0.98
x <- rep(0,1000)
x[1] <- rnorm(1,sd=omega)
for(tt in 2:1000) x[tt] <- rnorm(1,sd=sqrt(omega+alpha*x[tt-1]^2))
x <- round(x,6)
write.csv(x,file="arch1bsp.csv",row.names=F)

# logearnings.csv
set.seed(1234)
sigma.beta <- 0.02
rho <- 0.8
sigma.eta <- 0.15
NN <- 2000
TT <- 15
x <- matrix(NA,NN,TT)
beta.i <- rnorm(NN,sd=sigma.beta)
u <- rnorm(NN,sd=sigma.eta/sqrt(1-rho^2))
x[,1] <- beta.i+u
for(tt in 2:TT) {
  u <- rho*u+rnorm(NN,sd=sigma.eta)
  x[,tt] <- beta.i*tt+u
}
write.csv(round(x,5),file="logearnings.csv",row.names=F)

# acd1bsp.csv
setwd("c:/temp")
set.seed(1234)
omega <- 0.006
alpha <- 0.94
TT <- 2000
x <- rep(NA,TT)
x[1] <- rexp(1)*omega/(1-alpha)
for(tt in 2:TT) {
  psi <- omega+alpha*x[tt-1]
  x[tt] <- rexp(1)*psi
}
write.csv(round(x,6),file="acd1bsp.csv",row.names=F)

# oupath.csv
library(sde)
set.seed(1234)
mu <- 9
lambda <- 1.3
sigma <- 3
x <- sde.sim(t0=0,T=100,X0=mu, N=100,theta=c(lambda*mu,lambda,sigma),model="OU")
write.csv(round(x[-1],5),file="oupath.csv",row.names=F)

# timeaggr.csv
library(sde)
set.seed(123)
mu <- 0.00025
sigma <- 0.015
TT <- 360
x <- GBM(x=1,r=mu,sigma=sigma,T=TT,N=TT*100)
X <- matrix(x[-1],TT,100,byrow=T)
Y <- apply(X,1,sum)*0.01
arima(Y,order=c(1,0,1))
write.csv(round(Y,5),file="timeaggr.csv",row.names=F)

# rwnoise-csv
library(dlm)
set.seed(1234)
VV <- 9
WW <- 1
TT <- 1000
y <- th <- rep(NA,TT)
theta <- 5
for(tt in 1:TT) {
  theta <- theta+rnorm(1,sd=sqrt(WW))
  th[tt] <- theta
  y[tt] <- theta+rnorm(1,sd=sqrt(VV))
}
plot(y,t="l")
write.csv(round(y,5),file="rwnoise.csv",row.names=F)

# ttestboot.csv
set.seed(12345)
x <- c(1,9,2,8,3,7,4,6,5)
u <- rnorm(n=9,sd=2)
y <- round(1+u,6)
a <- data.frame(x=x,y=y)
write.csv(a,file="ttestboot.csv",row.names=F)

# olsgmm.csv
# Simulate data
TT <- 100; alph <- 1; bet1 <- 0.5; bet2 <- 1.2; sigm <- 0.5
u <- rnorm(TT,mean=0,sd=sigm)
x1 <- runif(TT,min=0,max=40)
x2 <- rnorm(TT,mean=10,sd=10)
y <- alph+bet1*x1 +bet2*x2 + u
z <- data.frame(y=round(y,2),x1=round(x1,2),x2=round(x2,2))
write.csv(z,file="olsgmm.csv",row.names=F)

# ar1bsp.csv
mu <- 5
rho <- 0.9
sigma2eps <- 1
TT <- 100
x <- rep(NA,TT)
x[1] <- rnorm(1,sd=sqrt(sigma2eps/(1-rho^2)))
for(i in 2:TT){
  x[i] <- mu+rho*(x[i-1]-mu)+rnorm(1,sd=sqrt(sigma2eps))
}
x <- data.frame(x=round(x,3))
write.csv(x,file="ar1bsp.csv",row.names=F)
