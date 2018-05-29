#' @title Root Mean Square Error (RMSE)
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description Estimates the Root Mean Square Error of two univariate signals Y (imputed values) and X (true values).
#' @details
#' This function returns the value of RMSE of two vectors corresponding to univariate signals.
#' A lower RMSE (\eqn{RMSE \in [0, \inf]}) value indicates a better performance method for the imputation task.
#' Both vectors Y and X must be of equal length, on the contrary an error will be displayed.
#' In both input vectors, eventual NA will be exluded with a warning diplayed.
#' @param Y vector of imputed values
#' @param X vector of true values
#' @examples
#' data(dataDTWBI)
#' X <- dataDTWBI[, 1] ; Y <- dataDTWBI[, 2]
#' compute.rmse(Y,X)

compute.rmse<-function(Y,X){
  
  if(length(Y)!=length(X)){stop("Input vectors are of different length !!!")}
  
  lengthNAX <- sum(is.na(X)) # Number of NA values
  if(lengthNAX > 0){warning(paste("Vector of true values contains ", lengthNAX, " NA !!! NA excluded", sep = ""))}
  lengthNAY <- sum(is.na(Y)) # Number of NA values
  if(lengthNAY > 0){warning(paste("Vector of imputed values contains ", lengthNAY, " NA !!! NA excluded", sep = ""))}
  n <- length(X)-max(lengthNAX, lengthNAY)
  
  out=sqrt(sum((Y-X)^2, na.rm = T)/n)
  
  return(out)
  
}