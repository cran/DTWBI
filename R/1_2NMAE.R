#' @title Normalized Mean Absolute Error (NMAE)
#' @author Camille Dezecache, Hong T. T. Phan, Emilie Poisson-Caillault
#' @description Estimates the Normalized Mean Absolute Error of two univariate signals Y (imputed values) and X (true values).
#' @details
#' This function returns the value of NMAE of two vectors corresponding to univariate signals.
#' A lower NMAE (\eqn{NMAE \in [0, \inf]}) value indicates a better performance method for the imputation task.
#' Both vectors Y and X must be of equal length, on the contrary an error will be displayed.
#' In both input vectors, eventual NA will be exluded with a warning diplayed.
#' @param Y vector of imputed values
#' @param X vector of true values
#' @examples
#' data(dataDTWBI)
#' X <- dataDTWBI[, 1] ; Y <- dataDTWBI[, 2]
#' compute.nmae(Y,X)
#' 
#' # If true values is a constant vector, NMAE = Inf.
#' # A warning is displayed and MAE is estimated instead of NMAE,
#' # unless true and imputed values are equal. In this case,
#' # by definition, NMAE = 0.
#' X <- rep(0, 10)
#' Y <- runif(10)
#' compute.nmae(Y,X) # MAE computed
#' Y <- X
#' compute.nmae(Y,X) # By definition, NMAE = 0

compute.nmae<-function(Y,X){
  
  if(length(Y)!=length(X)){stop("Input vectors are of different length !!!")}
  
  lengthNAX <- sum(is.na(X)) # Number of NA values
  if(lengthNAX > 0){warning(paste("Vector of true values contains ", lengthNAX, " NA !!! NA excluded", sep = ""))}
  lengthNAY <- sum(is.na(Y)) # Number of NA values
  if(lengthNAY > 0){warning(paste("Vector of imputed values contains ", lengthNAY, " NA !!! NA excluded", sep = ""))}
  n <- length(X)-max(lengthNAX, lengthNAY)

  if((max(X, na.rm = T)-min(X, na.rm = T))==0){
    if((max(Y, na.rm = T)-min(Y, na.rm = T))==(max(X, na.rm = T)-min(X, na.rm = T))){
      warning("Vectors of true and imputed values are constant and equal !!! By definition NMAE=0")
      out <- 0
      return(out)
    }
    else{
    warning("Vector of true values is constant !!! MAE was computed instead of NMAE !!!")
    numerator= sum(abs(Y-X), na.rm = T)
    out=numerator/n
    return(out)
    }
  }
  else{
  numerator= sum(abs(Y-X)/(max(X, na.rm = T)-min(X, na.rm = T)), na.rm = T)
  out=numerator/n
  return(out)
  }
}