# ===========================================================================
#
#      Stochastic Volatility
#
# ===========================================================================
rm (list = ls(all=TRUE))
graphics.off()


############
# Part (a) #
############
#Simulation estimation procedures circumvent this misspecification problem by 
#simulating the continuous model, and constructing a discrete time series of 
#simulated data which are calibrated with the observed discrete data thorough 
#a set of moment conditions based on a chosen auxiliary model. In an indirect
#inference approach these moment conditions are the parameters of the auxiliary
#model. The goal is to match these as close as possible.

############
# Part (b) #
############
# Function to generate discrete data from continuous time stochastic volatility model
svdata <- function(theta,y0,sig0,U,dt){
  n  <- dim(U)[1]            #     Length of the simulated series
  TT <- n*dt                 #     Discrete Sample size
  mu <- theta[1]; alpha <- theta[2]; kappa <- theta[3]; gama <- theta[4]
  expterm <- exp(-alpha*dt)
  bsigm <- expterm
  asigm <- kappa*(1-expterm)
  csigm <- gama*sqrt(dt)*sqrt((1-expterm*expterm)/(2*alpha))
  lnsig2 <- filter(asigm + csigm*U[,1],bsigm,method="recursive",init=log(sig0^2))
  sig2 <- exp(lnsig2)
  ay <- (mu-1/2*sig2)*dt
  by <- 1
  cy <- sqrt(sig2*dt)
  lny <- filter(ay + cy*U[,2],by,method="recursive",init=log(y0))
  lnyt <- rep(0, TT)
  for (t in 1:TT) {
    lnyt[t] <- lny[t*(1/dt)]
  }       
  return(lnyt)
}

TT         <- 100              #     Sample size                         
dt         <- 0.1              #     Continuous time step interval
n          <- TT/dt            #     Length of the simulated series
mu         <- 1                #     Mean parameter of geometric Brownian process
alpha      <- 0.9              #     Sensitivity parameter of mean reversion of Ornstein-Uhlenbeck process 
kappa      <- 0                #     Mean of Ornstein-Uhlenbeck process
gama       <- 1                #     sd of Ornstein-Ulenbeck process shock
y0         <- mu               #     Initial value of geometric Brownian process   
sig0       <- sqrt(exp(kappa)) #     Initial value of OU process

theta_true <- c(mu, alpha, kappa, gama)
U <- cbind(rnorm(n),rnorm(n))
lnytrue <- svdata(theta_true,y0,sig0,U,dt)
plot(lnytrue)


############
# Part (c) #
############
# Estimators:
# muhat <- mean(lny[2:TT] - lny[1:TT-1])
# rt <- lny[2:TT] - lny[1:TT-1] - muhat
# lnrt2 <- log(rt^2)
# ols <- lm(lnrt2[2:(TT-1)]~lnrt2[1:(TT-2)])
# deltahat0 <- ols$coefficients[1]
# deltahat1 <- ols$coefficients[2]
# deltahat2 <- sd(ols$residuals)
# thetahat <- rbind(muhat,deltahat0,deltahat1,deltahat2)
# delta0 = alpha*kappa
# delta1 = 1-alpha
# epsilon = ln(u1[2:TT]^2)-(1-alpha)*ln(u1[1:TT-1]^2)+gama*u2[2:TT]
# delta2 is sd(epsilon)

############
# Part (d) #
############
svind <- function(lnytrue,H,W,startval){
  # Estimate the auxiliary model using true data
  muhat <- mean(lnytrue[2:TT] - lnytrue[1:TT-1])
  rt <- lnytrue[2:TT] - lnytrue[1:TT-1] - muhat
  lnrt2 <- log(rt^2)
  ols <- lm(lnrt2[2:(TT-1)]~lnrt2[1:(TT-2)])
  deltahat0 <- ols$coefficients[1]
  deltahat1 <- ols$coefficients[2]
  deltahat2 <- sd(ols$residuals)
  thetahat <- rbind(muhat,deltahat0,deltahat1,deltahat2)
  
  # Generate std normal errors to simulate the true model      
  u1 <- matrix(rnorm(n*H), ncol=H,nrow=n)
  u2 <- matrix(rnorm(n*H), ncol=H,nrow=n)
  fobj <- function(theta) {
    thetahatsim <- matrix(NA,4,H)
    for (h in 1:H){
      #simulated data depending on theta
      lnysim <- svdata(theta,y0,sig0,cbind(u1[,h],u2[,h]),dt)
      muhatsim <- mean(lnysim[2:TT] - lnysim[1:TT-1])
      rtsim <- lnysim[2:TT] - lnysim[1:TT-1] - muhatsim
      lnrt2sim <- log(rtsim^2)
      olssim <- lm(lnrt2sim[2:(TT-1)]~lnrt2sim[1:(TT-2)])
      deltahat0sim <- olssim$coefficients[1]
      deltahat1sim <- olssim$coefficients[2]
      deltahat2sim <- sd(olssim$residuals)
      thetahatsim[,h] <- rbind(muhatsim,deltahat0sim,deltahat1sim,deltahat2sim)
    }
    thetatilde <- rowMeans(thetahatsim)
    Q <- t(thetahat - thetatilde) %*% W %*% (thetahat - thetatilde)
    return(Q)
  }  
  # Estimate the indirect model using simulated data   
  estResults <- optim(startval, fobj, method="L-BFGS-B", lower=c(0,0.1,-3,0.1), upper=c(3,3,3,3))
  return(estResults$par)
}

H          <- 10             #     Simulation runs
W          <- diag(4)        #     Optimal weight matrix in identified model
Z          <- svind(lnytrue,H,W,theta_true)
cat('\nTrue population parameter           = ', theta_true)
cat('\nEstimate                            = ', Z)

# Part e
R          <- 100            #     Draws to generate sampling disribution
Z <- matrix(NA,R,4)
pb <- txtProgressBar(min=0,max=R,style=3)
for (r in 1:R){
  setTxtProgressBar(pb,r)
  lnytrue <- svdata(theta_true,y0,sig0,cbind(rnorm(n),rnorm(n)),dt)
  Z[r,] <- svind(lnytrue,H,W,theta_true)
}
close(pb)

# Generate statistics on the sampling distribution    
cat('\n')
cat('\nNumber of replications              = ', R)
cat('\nSample size                         = ', TT)
cat('\n')
cat('\nTrue population parameter           = ', theta_true)
cat('\nMean of estimates                   = ', colMeans(Z))
cat('\nStandard deviation of estimates     = ', apply(Z,2,sd))
cat('\nRMSE of Theta                       = ', sqrt(colMeans(sweep(Z,2,theta_true)^2)))
hist(Z[,1],prob=T,main = "Histogram of mu")
hist(Z[,2],prob=T,main = "Histogram of alpha")
hist(Z[,3],prob=T,main = "Histogram of kappa")
hist(Z[,4],prob=T,main = "Histogram of gamma")


############
# Part (e) #
############
# Since this is a just identified model, the weight matrix does not matter.
#For over-identified models, identity matrix as weight matrix means that indirect
#inference estimator is still consistent but not asymptotically efficient.