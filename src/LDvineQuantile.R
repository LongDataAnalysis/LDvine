#' @title Conditional quantile prediction function for 5-dimesional D-vine model 
#'
#' @description This function calculates Conditional quantile prediction for 5-dimensional D-vine model.
#'
#' @param family.set list of families for a D-Vine
#' @param p specified quantile level
#' @param par.set list of fitted copula parameter values
#' @param u dataframe that contains the data information of the measurements
#' @return estimated quantiles on the copula scale 
#' @examples
#'
#' # Define vine components as a list
#'
#' family.set =  list(c(4,4,4,4),c(14,1,13),c(5,1),c(5))
#' par.set = list(c(1.2420, 1.2524, 1.2865, 1.4522), c(1.2287, 0.2011, 0.3443), c(0.9139, 0.2964), c(0.2045))
#' p=0.5
#' data=data.frame(0.3020696,0.9416951, 0.9032743, 0.3863236)
#' u_quant = DvineQuant (p=0.5, u=data, family.set = family.set , par.set = par.set)
#' 
#' 
#' @export

############################################################################ 
# Function to calculate Conditional Quantile 
############################################################################
 DvineQuant = function(p, u, family.set, par.set){

  # obtain the true conditional marginals (Tree 2)
  vv1.2 = VineCopula::BiCopHfunc2(u1=u[,1],u2=u[,2],family = family.set[[1]][1],par=par.set[[1]][1])
  vv3.2 = VineCopula::BiCopHfunc2(u1=u[,3],u2=u[,2],family = family.set[[1]][2],par=par.set[[1]][2])
  
  vv2.3 = VineCopula::BiCopHfunc2(u1=u[,2],u2=u[,3],family = family.set[[1]][2],par=par.set[[1]][2])
  vv4.3 = VineCopula::BiCopHfunc2(u1=u[,4],u2=u[,3],family = family.set[[1]][3],par=par.set[[1]][3])
  
  # obtain the true conditional marginals (Tree 3)
  vvv1.23 = VineCopula::BiCopHfunc2(vv1.2,  vv3.2,  family = family.set[[2]][1], par=par.set[[2]][1])
  vvv4.23 = VineCopula::BiCopHfunc2(vv4.3,  vv2.3,  family = family.set[[2]][2], par=par.set[[2]][2])
  
  # obtain the true conditional marginals (Tree 4)
  vvvv1.234 = VineCopula::BiCopHfunc2(vvv1.23,  vvv4.23,  family = family.set[[3]][1], par=par.set[[3]][1])
  

  #############################################################################
  # Conditional quantile function in copula scale
  #############################################################################
  
  # first h-inverse function
  
  hh5.1234 <- VineCopula::BiCopHinv2(p, vvvv1.234,family=family.set[[4]][1], par=par.set[[4]][1])
  
  # second h-inverse function
  vvv2.34 = VineCopula::BiCopHfunc2(vv2.3,  vv4.3,  family = family.set[[2]][2], par=par.set[[2]][2])
  
  hh5.234 <- VineCopula::BiCopHinv2(hh5.1234, vvv2.34,family=family.set[[3]][2], par=par.set[[3]][2])
  
  # third h-inverse function
  vv3.4 = VineCopula::BiCopHfunc2(u1=u[,3],u2=u[,4],family = family.set[[1]][3], par=par.set[[1]][3])
  
  hh5.34 <- VineCopula::BiCopHinv2(hh5.234, vv3.4,family=family.set[[2]][3], par=par.set[[2]][3])
  
  # final h-inverse function
  hh5.4 <- VineCopula::BiCopHinv2(hh5.34, u[,4],family = family.set[[1]][4], par=par.set[[1]][4])
  
  res<- hh5.4
  return(res)
  
}

 