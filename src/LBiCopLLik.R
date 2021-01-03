#' @title Log-likelihood calculation in bivariate copulas for longitudinal data
#'
#' @description This function calculates log-likelihood function for decay parameters in bivariate copulas.
#'
#' @param u1,u2 numeric vectors of equal length with values in $[0,1]$.
#' @param x numeric value of time interval between two visits.
#' @param family numeric value of the parametric copula family as coded in VineCopula package.
#' @param beta vector of numeric values of the decay parameters.
#' @param repar specifies whether decay parameterization is done in terms of Pearson's rho or Kendall's tau.
#' @return returns the log-likelihood values.
#' @example 
#' LBiCopLLik(u1=0.2, u2=0.3, x=0.15, family=3, beta=c(0.3, 0.2), repar="pearson")
#' LBiCopLLik(u1=runif(1,0,1), u2=runif(1,0,1), x=runif(1,0,1), family=3, beta=c(0.3, 0.2), repar="pearson")
#' @export
LBiCopLLik  <- function(u1, u2, x, family, beta, repar=c("kendall", "pearson", "spearman")){

  npar <- length(beta)
  
  # Reparameterization:
  repar <- match.arg(repar)

  if(repar == "kendall"){
    if(npar==1){tau <- exp(-beta*x)}
    if(npar==2){tau <- exp(-beta[1] - beta[2]*x) }
    tau <- ifelse(tau>0.96, 0.96, tau)
    par <- VineCopula::BiCopTau2Par(family=family, tau=tau)
  }

  if(repar == "pearson"){
    if(npar==1){rho <- exp(-beta*x)}
    if(npar==2){rho <- exp(-beta[1] - beta[2]*x) }
    tau <- 2/pi*(asin(rho))
    tau <- ifelse(tau>0.96, 0.96, tau)
    tau <- ifelse(tau==0, 0.001, tau)
    par <- VineCopula::BiCopTau2Par(family=family, tau=tau)
    par <- ifelse(par==0, 0.001, par) # to control conversion problem.
  }

  if(repar == "spearman"){
    #srho <- exp(-beta*x)
    if(npar==1){srho <- exp(-beta*x)}
    if(npar==2){srho <- exp(-beta[1] - beta[2]*x) }
    rho <- 2*sin(srho*pi/6)
    tau <- 2/pi*(asin(rho))
    tau <- ifelse(tau>0.96, 0.96, tau)
    par <- VineCopula::BiCopTau2Par(family=family, tau=tau)
  }

  dens <- VineCopula::BiCopPDF(u1,u2, family=family, par=par)
  ldens <- log(dens)

 return(loglik=sum(ldens))

}



