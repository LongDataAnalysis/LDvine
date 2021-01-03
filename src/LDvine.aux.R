############################################################################
## AUXILLARY FUNCTIONS FOR SIMULATION OF TIME-DEPENDENT LONGITUDINAL DATA ##     
############################################################################

############################################################################ 
# Function to calculate rho matrix
############################################################################
# Tau matrix for a subject with visits given by the Time vector, where the dependence decays exponentially with parameter beta

LDVine.Tau <- function(beta0,beta, Time){
  d <- length(Time)
  Tau <- matrix(0,nrow=d,ncol=d)
  rho <- matrix(0,nrow=d,ncol=d)
  
  # If beta is a scalar it is converted to a matrix:
  if(length(beta0)==1){ 
    alpha0 <- sapply(1:d, function(k) c(rep(0,k), rep(beta0,(d-k)) )) 
  }
  
  if(length(beta)==1){ 
    alpha <- sapply(1:d, function(k) c(rep(0,k), rep(beta,(d-k)) )) 
  }
  if(!is.matrix(beta)) stop("Provide a ", d, " x ", d, " lower triangular matrix containing beta parameters.")
  
  
  for(k in d:2){
    rho[k, 1:(k-1)] <- exp(-beta0[k,1:(k-1)] -beta[k,1:(k-1)] * diff(Time, lag = (d-k+1)))
    Tau[k, 1:(k-1)] <- 2/pi*(asin(rho[k, 1:(k-1)])) 
  }
  return(Tau)
}


############################################################################ 
# Function to convert rho to TAU
############################################################################
# Convert Tau matrix to Par matrix given an RVM matrix of copula families

LDVine.rho2tau <- function(beta, x){
  tau = 2/pi*asin(exp(-beta[1]-beta[2]*x))
  return(tau)
}  

############################################################################ 
# Function to convert TAU to PAR
############################################################################
# Convert Tau matrix to Par matrix given an RVM matrix of copula families

LDVine.Par <- function(Family, Tau){
  Par = VineCopula::BiCopTau2Par(Family, Tau)
  return(Par)
} 

############################################################################ 
# Function to generate TIME variable
############################################################################
# Generation of visit times for one subject

Longt.Sim<-function(d){
  Time<-sapply(1:d, function(k){runif(1, (k-1)/d, k/d)})
  return(Time)
}


############################################################################ 
# Function to Simulate time-heterogeneous longitudinal data
############################################################################

#' @title Simulation of time-heterogeneous longitudinal data.
#'
#' @description This function simulates time-heterogeneous longitudinal copula data from D-vine copula.
#'
#' @param n sample size
#' @param d dimension
#' @param Family list of families for a D-Vine
#' @param par list of parameter functions for a D-Vine
#' @param par2 list of second parameters for two dimensional copula families 
#' @return A data matrix having non-iid distribution from a vine copula. 
#' @examples
#'
#' # Define RVM object
#'
#' library(VineCopula)
#' Family = matrix(c(0,5,3,4,0,0,3,4,0,0,0,4,0,0,0,0),4,4)
#' beta0 = matrix(c(0, 0.8, 0.7, 0.6,0,0,0.7,0.6,0,0,0,0.6,0,0,0,0),4,4)
#' beta = matrix(c(0, 0.8, 0.7, 0.6,0,0,0.7,0.6,0,0,0,0.6,0,0,0,0),4,4)
#' Time = t(sapply(1:100, function(i) Longt.Sim(4)))
#' Tau.Array = sapply(1:100, function(i) LDVine.Tau(beta0, beta, Time[i,]), simplify="array")
#' Par.Array = sapply(1:100, function(i) LDVine.Par(Family, Tau = Tau.Array[,,i]),  simplify="array")
#' par2  =  matrix(rep(0,16),4,4)
#' data = data.frame(t(sapply(1:100, function(i)  LDVine.Sim(d=4, Family, Par.Array[,,i]))))
#' 
#' @export

# Data generation for a subject given an RVM matrix for copula families and Par matrix based on the times of d visits
LDVine.Sim <- function(d, Family, Par, Par2 = matrix(0,nrow=d,ncol=d)){ 
  
  if(!is.matrix(Family)) stop("Provide a lower triangular matrix containing copula families.")
  if(!is.matrix(Par)) stop("Provide a lower triangular matrix containing copula parameter values.")
  if(nrow(Par) != d | nrow(Family) != d) stop("Tau and Family should be d x d matrices.")
  
  elements <- 1:d
  Mat <- diag(1:d)
  for (i in 2:d){
    Mat[d+2-i, 1:(d+1-i)] <- elements[i:d]
  }
  
  RVM <- VineCopula::RVineMatrix(Mat, Family, Par, Par2)
  U <- VineCopula::RVineSim(N=1, RVM = RVM)
  return(U)
}
