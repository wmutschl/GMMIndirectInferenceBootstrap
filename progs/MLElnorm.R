# Schätzung der zensierten Lognormalverteilung

# Ziehen der Stichprobe
x <- rlnorm(1000,meanlog=1,sdlog=0.5)
bbg <- 5
x[x>bbg] <- bbg # Zensierung an der Stelle cc

# Aufstellen der Loglikelihoodfunktion
logliklnorm <- function(param,cc,daten)
  {
  mu <- param[1]
  sigma2 <- param[2]
  c1 <- sum(log(dlnorm(daten[daten<cc],meanlog=mu,sdlog=sqrt(sigma2))/plnorm(cc,meanlog=mu,sdlog=sqrt(sigma2))))
  c2 <- sum(log(1-plnorm(cc,meanlog=mu,sdlog=sqrt(sigma2))))
  return(c1+c2)
  }

# Numerische Optimierung
mu0 <- 0.5 # Startwerte
sigma20 <- 0.3
opt1 <- optim(c(mu0,sigma20),logliklnorm,control=list(fnscale=-1),hessian=TRUE,cc=bbg,daten=x)

muhat <- opt1$par[1]
sigma2hat <- opt1$par[2]
covhat <- -solve(opt1$hessian)
