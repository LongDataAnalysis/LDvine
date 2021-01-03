#' @title Estimation of decay parameter in bivariate copulas
#'
#' @description This function performs maximum likelihood estimation of the decay parameter for bivariate copulas.
#'
#' @param u1,u2 numeric vectors of equal length with values in $[0,1]$.
#' @param x numeric value of time interval between two visits.
#' @param family numeric vector specifying parametric copula families.
#' 1 = Gaussian copula
#'
#' 2 = Student t copula (t-copula)
#'
#' 3 = Clayton copula
#'
#' 4 = Gumbel copula
#'
#' 5 = Frank copula
#'
#' 6 = Joe copula
#'
#' 7 = BB1 copula
#'
#' 8 = BB6 copula
#'
#' 9 = BB7 copula
#'
#' 10 = BB8 copula
#'
#' 13 = rotated Clayton copula (180 degrees; “survival Clayton”)
#'
#' 14 = rotated Gumbel copula (180 degrees; “survival Gumbel”)
#'
#' 16 = rotated Joe copula (180 degrees; “survival Joe”)
#'
#' 17 = rotated BB1 copula (180 degrees; “survival BB1”)
#'
#' 18 = rotated BB6 copula (180 degrees; “survival BB6”)
#'
#' 19 = rotated BB7 copula (180 degrees; “survival BB7”)
#'
#' 20 = rotated BB8 copula (180 degrees; “survival BB8”)
#'
#' 23 = rotated Clayton copula (90 degrees)
#'
#' 24 = rotated Gumbel copula (90 degrees)
#'
#' 26 = rotated Joe copula (90 degrees)
#'
#' 27 = rotated BB1 copula (90 degrees)
#'
#' 28 = rotated BB6 copula (90 degrees)
#'
#' 29 = rotated BB7 copula (90 degrees)
#'
#' 30 = rotated BB8 copula (90 degrees)
#'
#' 33 = rotated Clayton copula (270 degrees)
#'
#' 34 = rotated Gumbel copula (270 degrees)
#'
#' 36 = rotated Joe copula (270 degrees)
#'
#' 37 = rotated BB1 copula (270 degrees)
#'
#' 38 = rotated BB6 copula (270 degrees)
#'
#' 39 = rotated BB7 copula (270 degrees)
#'
#' 40 = rotated BB8 copula (270 degrees)
#'
#' 104 = Tawn type 1 copula
#'
#' 114 = rotated Tawn type 1 copula (180 degrees)
#'
#' 124 = rotated Tawn type 1 copula (90 degrees)
#'
#' 134 = rotated Tawn type 1 copula (270 degrees)
#'
#' 204 = Tawn type 2 copula
#'
#' 214 = rotated Tawn type 2 copula (180 degrees)
#'
#' 224 = rotated Tawn type 2 copula (90 degrees)
#'
#' 234 = rotated Tawn type 2 copula (270 degrees)
#' @param repar specifies whether decay parameterization is done in terms of Pearson's rho or Kendall's tau.
#' @param start specifies the initial value of parameters
#' @return returns the decay parameter estimates with selected family based on minimum AIC.
#' @examples
#' LBiCopFit(u1=0.2, u2=0.4, x=1.25, family=c(3,4), repar="pearson", start=c(0.1,0.1))
#' LBiCopFit(u1=runif(1,0,1), u2=runif(1,0,1), x=runif(1,0,1), family=c(3,4,5), repar="pearson", start=c(0.1,0.1))
#' @export

LBiCopFit <- function(u1, u2, x, family, repar=c("kendall", "pearson", "spearman"), start){

  repar <- match.arg(repar)
  res <- vector("list", length(family))
  aic <- rep(NA, length(family))
  
  # subsetting the data to extract the available observations (in case of unbalanced data)
  U.data<-data.frame(u1,u2)
  index<-which(complete.cases(U.data))
  data<-U.data[index,]
  for (k in seq_along(family)){
    obj.fnc <- function(beta){-LBiCopLLik(data[,1], data[,2], x[index], family=family[k], beta, repar)}
    res[[k]]<- nlminb(start = start , objective = obj.fnc, lower = 0.01)
    aic[k] <- -2*(-res[[k]]$objective) + 2*2 
  }
  # Select family based on minimum AIC:
  ch <- which.min(aic)
  family <- family[ch]
  fit <- unlist(res[[ch]])
  
  #if (length(start)==1){return(list(beta = as.numeric(fit[1]), loglik = -as.numeric(fit[2]), family=family))}
  if (length(start)==2){return(list(beta0 = as.numeric(fit[1]), beta1 = as.numeric(fit[2]), loglik = -as.numeric(fit[3]), family=family))}
}


