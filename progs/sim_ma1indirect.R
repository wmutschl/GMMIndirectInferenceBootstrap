#============================================================================
#
#      Monte Carlo analysis to investigate the sampling properties
#      of the indirect estimator of a first order MA model.
#
#      Gourieroux et. al. (1993) J of Appl. Eco.
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()
set.seed(123, kind="Mersenne-Twister")

#
#------------------------- Helper Functions----------------------------------
#

# load required functions - recserar, trimr
source("EMTSUtil.R")

#-------------------------------------------------------------------------#
# Objective function to compute the indirect estimator 
#-------------------------------------------------------------------------#
q <- function(b,e,lag,bhat) {
  # Simulate data making sure that b[1] is in the unit circle 
  ys     <- trimr(e,1,0) - tanh(b)*trimr(e,0,1)        
  ys     <- ys - mean(ys)
  
  if (lag == 1)
    bhats <- lm(trimr(ys,1,0) ~ trimr(ys,0,1) - 1)$coef
  else if (lag == 2)
    bhats = lm(trimr(ys,2,0) ~ cbind(trimr(ys,1,1),trimr(ys,0,2)) - 1 )$coef
  else if (lag == 3)
    bhats <- lm(trimr(ys,3,0) ~ cbind(trimr(ys,2,1),trimr(ys,1,2),trimr(ys,0,3)) - 1)$coef
  val <- t( (bhat - bhats) ) %*% (bhat - bhats)
  return(val) 
}

# 
#------------------- Indirect Estimation of the AR(1) Model -----------------
#

sim_ma1indirect <- function() 
{
    # Parameters of MA(1) process
  t     <- 250                                
  theta <- 0.5                    
  lag   <- 3             

  # Simulation settings  
  nreps <- 1000
  b  <- rep(0, nreps)

  # Main DO LOOP to generate sampling disribution
  pb <- txtProgressBar(min=0, max=nreps, style=3)
  for (j in seq(nreps)) {
    # Generate the actual data for the MA(1) process       
    u  <- rnorm(t)
    y  <- trimr(u,1,0) - theta*trimr(u,0,1)
    
    # Estimate the the auxiliary model using actual data  
    y   <- y - mean(y)
    if (lag == 1)
      bhat <- lm(trimr(y,1,0) ~ trimr(y,0,1) - 1)$coef
    else if (lag == 2)
      bhat = lm(trimr(y,2,0) ~ cbind(trimr(y,1,1),trimr(y,0,2)) - 1 )$coef
    else if (lag == 3)
      bhat <- lm(trimr(y,3,0) ~ cbind(trimr(y,2,1),trimr(y,1,2),trimr(y,0,3)) - 1)$coef
    
    # Compute indirect estimator (could use a line search algorithm)
    e    <- rnorm(t)
    estResults <- optim(theta, q, e=e, lag=lag, bhat=bhat, method="BFGS")
    b[j] <- estResults$par    
    setTxtProgressBar(pb, j)    
  }
  close(pb)

  # Generate statistics on the sampling distribution
  b <- tanh(b)
  m     <- mean(b)
  stdev <- sd(b)
  rmse  <- sqrt(mean(b-theta)^2)

  cat('\n')
  cat('\nNumber of replications              =  ',nreps)
  cat('\nSample size                         =  ', t, '\n')
  print(cbind("True"=theta, "Mean"=m, "Std err"=stdev, "RMSE"=rmse))
}
