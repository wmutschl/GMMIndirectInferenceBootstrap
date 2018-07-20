# ===========================================================================
#
#      Finite sample properties for the Gamma model
#
# ===========================================================================

rm (list = ls(all=TRUE))
graphics.off()
library(gmm)

#-------------------------------------------------------------------------
# Negative log-likelihood function for gamma distribution with beta = 1
#-------------------------------------------------------------------------
neglog <- function(a,y) {
  b <- 1
  n <- length(y)
  logl <-  (a-1)*sum(log(y)) - sum(y)/b -n*a*log(b) - n*lgamma(a)
  return(-logl)
}

#-------------------------------------------------------------------------
# GMM moment function for gamma distribution based on first moment
#-------------------------------------------------------------------------
g1 <- function(a,y) {
  b<-1
  m <- (y-a*b)
  return(cbind(m)) 
}

#-------------------------------------------------------------------------
# GMM moment function for gamma distribution based on first two moments
#-------------------------------------------------------------------------
g2 <- function(a,y) {
  b <- 1
  m1 <- (y-a*b)
  m2 <- (y^2 - a^2 - a)
  m <- cbind(m1,m2)
  return(m) 
}

#-------------------------------------------------------------------------
#  Monte Carlo experiment for different parameter inputs
#-------------------------------------------------------------------------

simexp <- function(TT,reps,dist,a0,H0) {
  #TT:sample size, reps:no of repititions, dist:gam or exp, a0: true parameter values

  cols <- length(a0) # how many different true values
  
  # Initialize output
  zeros <- array(0, c(reps, cols))
  # Point estimate
  aML  <- zeros; aGMM1   <- zeros; aGMM2   <- zeros
  # Bias
  biasML  <- zeros; biasGMM1   <- zeros; biasGMM2   <- zeros
  # Standard error
  sdML <- zeros; sdGMM1  <- zeros; sdGMM2  <- zeros
  # t statistic
  tML  <- zeros; tGMM1   <- zeros; tGMM2   <- zeros

  for (j in 1:cols) {
    print(paste("Parameter",paste(j,cols,sep=" out of ")))
    pb <- txtProgressBar(min=0,max=reps,style=3)
    for (r in 1:reps) {
      setTxtProgressBar(pb,r)
      # Generate data
      if (dist == "gam") {y <- rgamma(TT, shape=a0[j],scale=1)}
      if (dist == "exp") {
        y <- rexp(TT,rate=1/a0[j])
        #y <- -log(runif(TT))*a0[j]
        }
      
      # Maximum likelihood
      estResultsML <- optim(mean(y),neglog,y=y,method="Brent",hessian = T,lower=0,upper=10)
      aML[r, j] <- estResultsML$par
      HessML <- estResultsML$hessian
      biasML[r,j]  <- aML[r,j]-a0[j]
      sdML[r,j] <- sqrt(1/HessML)
      tML[r,j]  <- (aML[r,j]-H0)/sdML[r,j]
      
      # GMM using the first moment
      estResultsGMM1 <- gmm(g1,y,t0=mean(y),wmatrix="optimal",optfct = "nlminb",lower=0,upper=10)
      aGMM1[r, j] <- estResultsGMM1$coefficients
      biasGMM1[r,j]  <- aGMM1[r,j]-a0[j]
      sdGMM1[r,j] <- sqrt(estResultsGMM1$vcov)
      tGMM1[r,j]  <- (aGMM1[r,j]-H0)/sdGMM1[r,j]
      
      # GMM using two moments
      estResultsGMM2 <- gmm(g2,y,t0=mean(y),wmatrix="optimal",optfct = "nlminb",lower=0,upper=10)
      aGMM2[r, j] <- estResultsGMM2$coefficients
      biasGMM2[r,j]  <- aGMM2[r,j]-a0[j]
      sdGMM2[r,j] <- sqrt(estResultsGMM2$vcov)
      tGMM2[r,j]  <- (aGMM2[r,j]-H0)/sdGMM2[r,j]
    }
    close(pb)
  }
  # Save results
  part_a <- as.data.frame(cbind(colMeans(biasML),colMeans(sdML),colMeans(biasGMM1),colMeans(sdGMM1),colMeans(biasGMM2),colMeans(sdGMM2)),row.names=as.character(a0))
  names(part_a) <- c("Bias ML","Stderr ML","Bias GMM1","Stderr GMM1","Bias GMM2","Stderr GMM2")
  part_b <- as.data.frame(cbind(colMeans(abs(tML) >1.96),colMeans(abs(tGMM1)>1.96),colMeans(abs(tGMM2)>1.96)),row.names=as.character(a0))
  names(part_b) <- c("|t ML| > 1.96","|t GMM1| > 1.96","|t GMM2| > 1.96")
  return(list("part_a"=part_a,"part_b"=part_b))
}

############
# Part (a) #
############
reps <- 1000
TT <- c(50,100,200)
dist <- 'gam'
a0 <- c(1,3,5)
H0 <- 1
Za <- list()
for (i in 1:length(TT)){
  name <- paste("Sample size",TT[i])
  print(name)
  Za[[name]] <- simexp(TT[i],reps,dist,a0,H0)
}
(Za$`Sample size 50`$part_a)
(Za$`Sample size 100`$part_a)
(Za$`Sample size 200`$part_a)
# BIAS:
# The ML estimator shows very small positive bias that decreases as T increases
# The GMM1 estimator is exactly the sample mean which is unbiased, so any nonzero bias of GMM1 reflects simulation error
# GMM2 estimator shows a small negative bias that tends to increase if a0 increases, however it decreases in magnitue as T increases (consistency!)

# SD:
# The standard error of all three estimators increases as a0 increases
# The standard error of all three estimatros decreases as T increases
# The variances approximately halve as the sample size doubles, since estimators are based on the covmatrix T^{-1}Omega
# sd(ML)<sd(GMM2)<sd(GMM1)
# ML is based on full distribution, therefore it is most efficient, GMM2 is based on 2 moments therefore it is more efficient than GMM1



############
# Part (b) #
############
reps <- 1000
TT <- 200
dist <- 'gam'
a0 <- c(1.00,1.05,1.10,1.15,1.20,1.25,1.30)
H0 <- 1
name <- paste("Sample size 200")
print(name)
Zb <- simexp(TT,reps,dist,a0,H0)
(Zb$part_b)
# The size of the test is given for a0=1 and the power for a0>1
# ML and tGMM1 have finite sample sizes closest to the nominal size of 0.05
# tGMM2 is marginally oversized
# The power of ML is highest, since ML is more efficient. GMM2 has higher power than GMM1
# More efficient estimators translate to more powerful hypothesis tests


############
# Part (c) #
############
reps <- 1000
TT <- c(200,400)
dist <- 'exp'
a0 <- c(1,3,5)
H0 <- 1
Zc <- list()
for (i in 1:length(TT)){
  name <- paste("Sample size",TT[i])
  print(name)
  Zc[[name]] <- simexp(TT[i],reps,dist,a0,H0)
}
(Zc$`Sample size 200`$part_a)
(Zc$`Sample size 400`$part_a)

# Note that the GMM estimator based on the first moment is correctly specified, 
# whereas the ML and GMM estimator based on the first two moments are only specified 
# correctly for the special case of a0=1 since the exponential and gamma distribution then coincide.
# All three estimators have small bias for a0=1, since gamma distribution is correctly specified.
# ML and GMM2 become increasingly negatively biased for a0>1, biases do not disappear in large samples
# GMM1 is unbiased and consistent
# If distributional assumption is incorrect ML is no longer consistent
# Provided moments are specified correctly, GMM is more robust against misspecification
