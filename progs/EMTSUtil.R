#====================================================================
# A set of helper functions used with the EMTS Book.
#===================================================================

#--------------------------------------------------------------------
# Graphic device compatible function
# Opens a new window device based on the OS.
#--------------------------------------------------------------------
figure <- function(){
  if(Sys.info()["sysname"] == "Windows") {
    windows()  
  }else {
    x11()
  }  
} 

#--------------------------------------------------------------------
# Matrix inverse function
# Wrapper function for computing matrix inverse
#--------------------------------------------------------------------
inv <- function (A) {
  return(solve(A))
}

#--------------------------------------------------------------------
# Numeric hessian function
# Computes finite difference Hessian matrix
#--------------------------------------------------------------------

library("nlme")

numhess <- function ( f, ... ) {  
  # Use nlme package for the hessian
  H <- fdHess(..., fun=f)$Hessian
  return(H)
}

#--------------------------------------------------------------------
# Numeric gradient function
# Computes numerical gradient matrix at each observation
#--------------------------------------------------------------------
numgrad <- function( fun,x, ... ) {
  f0  <-  fun(x, ... )             # n by 1    
  n   <- length( f0 )
  k   <- length( x )
  fdf <- array(0, dim=c(n,k) )
 
  # Compute step size 
  eps <- .Machine$double.eps
  dx  <- sqrt( eps )*( abs( x ) + eps )
  xh  <- x + dx
  dx  <- xh - x
  ind <- dx < sqrt(eps)
  dx[ind] <- sqrt(eps)
  
  # Compute gradient
  xdx <- diag(dx) + x  
  for (i in seq(k)) {
    fdf[,i] <- fun(xdx[,i], ...)    
  }   
  G0 <- kronecker(matrix(1, 1, k), f0 )                       # n by k 
  G1 <- kronecker(matrix(1, n, 1), t(dx) )  
  G  <- ( fdf-G0 ) / G1
  return(G)  
}


#--------------------------------------------------------------------
# Computes a vector of autoregressive recursive series
#
#   Inputs: 
#             x  = a matrix of dimensions (n,k) (exogenous variables)
#             y0 = a matrix of dimensions (p,k) (starting values)
#             a  = a matrix of dimensions (p,k) (lag parameters)
#
#  Mimics the Gauss routine. 
#--------------------------------------------------------------------
recserar <- function(x,y0,a){
  rx <- nrow(x)
  cx <- ncol(x)
    
  ry <- nrow(y0)
  cy <- ncol(y0)
  
  ra <- nrow(a)
  ca <- ncol(a)
  nargin <- length(as.list(match.call())) - 1      
  if ((nargin != 3) 
        || (cx != cy) 
        || (cx != ca) 
        || (ry != ra)) {
     stop('Check function inputs')    
  }    
  y <- array(0, c(rx,cx))
  for (j in 1:ry){
    y[j,] <- y0[j,]    
  }
  for (j in (ry+1):rx){
    y[j,] <- x[j,]
    for (k in 1:ry) {
      y[j,] <- y[j,] + a[k,] * y[j-k,]
    }     
  }   
  return(y)
}

#--------------------------------------------------------------------
# Returns a matrix (or vector) stripped of the specified rows
#
#   Inputs: 
#             x  = input matrix (or vector) (n x k)
#             rb = first n1 rows to strip
#             re = last  n2 rows to strip
#
#--------------------------------------------------------------------
trimr <- function(x,rb,re) {
  x <- cbind(x)
  n <- nrow(x)
  if ((rb+re) >= n) {
      stop('Attempting to trim too much')
  }
  z <- x[(rb+1):(n-re),]
  return(z)  
}

#--------------------------------------------------------------------
# 
# Puts given vector (v) onto the main diagonal of a square matrix (x)
# Reproduces the GAUSS command diagrv
#--------------------------------------------------------------------

diagrv <- function(X, v) {
  diag(X) <- v
  return(X)
}

#----------------------------------------------------------------------------
# Reshape a matrix to agree with GAUSS reshape command
#----------------------------------------------------------------------------
reshapeg <- function(Y,r,c) {
  library("matlab")
   tmp <- reshape(t(Y),c,r)
   X   <- t(tmp)
   return(X)
}

#--------------------------------------------------------------------
#
#   Mimic GAUSS vech( ) function by stacking columns of x
#   on and below the diagonal
#
#--------------------------------------------------------------------
vech <- function(x) {
  rows <- nrow(x)
  cols <- ncol(x)  
  y <- c()
  for (i in seq(rows)) {
    y <- c(y, x[i,1:i])
  }
  y <- cbind(y)  
}

#--------------------------------------------------------------------
#
#   Program to mimic the GAUSS seqa command
#
#--------------------------------------------------------------------
seqa <- function(start,step,nvals) {
  y <- seq(from=start, by=step, length.out=nvals)
  return(y)  
}

#----------------------------------------------------------------------------
# A simple lag function for matrices
#-----------------------------------------------------------------------------
lag.matrix <- function(m, nLags) {
  nargin <- length(as.list(match.call())) - 1      
  if (nargin != 2) {        
     stop('Check function inputs')    
  }
  lagM <- c()
  for(i in seq(nLags)) {
    for(j in seq(ncol(m))) {
      tmp <- c(rep(NA, i), trimr(m[,j], 0, i))      
      lagM <- cbind(lagM, tmp)      
    }    
  }  
  return(lagM)
}

#----------------------------------------------------------------------------
## A: coefficient matrix
## tol: tolerance for checking for 0 pivot
## verbose: if TRUE, print intermediate steps
## fractions: try to express nonintegers as rational numbers
## Written by John Fox
#----------------------------------------------------------------------------
rref <- function(A, tol=sqrt(.Machine$double.eps),verbose=FALSE,
                 fractions=FALSE){
  
  if (fractions) {
    mass <- require(MASS)
    if (!mass) stop("fractions=TRUE needs MASS package")
  }
  if ((!is.matrix(A)) || (!is.numeric(A)))
    stop("argument must be a numeric matrix")
  n <- nrow(A)
  m <- ncol(A)
  for (i in 1:min(c(m, n))){
    col <- A[,i]
    col[1:n < i] <- 0
    # find maximum pivot in current column at or below current row
    which <- which.max(abs(col))
    pivot <- A[which, i]
    if (abs(pivot) <= tol) next     # check for 0 pivot
    if (which > i) A[c(i, which),] <- A[c(which, i),]  # exchange rows
    A[i,] <- A[i,]/pivot            # pivot
    row <- A[i,]
    A <- A - outer(A[,i], row)      # sweep
    A[i,] <- row                    # restore current row
    if (verbose)
      if (fractions) print(fractions(A))
    else print(round(A,round(abs(log(tol,10)))))
  }
  for (i in 1:n)
    if (max(abs(A[i,1:m])) <= tol)
      A[c(i,n),] <- A[c(n,i),] # 0 rows to bottom
  if (fractions) fractions (A)
  else round(A, round(abs(log(tol,10))))
}


    
