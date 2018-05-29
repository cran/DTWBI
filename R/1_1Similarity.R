#' @title Similarity
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description Estimates the percentage of similarity of two univariate signals Y (imputed values) and X (true values).
#' @details
#' This function returns the value of similarity of two vectors corresponding to univariate signals.
#' A higher similarity (\eqn{Similarity \in [0, 1]}) highlights a more accurate method for completing missing values in univariate datasets.
#' Both vectors Y and X must be of equal length, on the contrary an error will be displayed.
#' In both input vectors, eventual NA will be excluded with a warning diplayed.
#' @param Y vector of imputed values
#' @param X vector of true values
#' @examples
#' data(dataDTWBI)
#' X <- dataDTWBI[, 1] ; Y <- dataDTWBI[, 2]
#' compute.sim(Y,X)
#' 
#' # By definition, if true values is a constant vector
#' # and one or more imputed values are equal to the true values,
#' # similarity = 1.
#' X <- rep(2, 10)
#' Y <- X
#' compute.sim(Y,X)

compute.sim<-function(Y,X){
  
  if(length(Y)!=length(X)){stop("Input vectors are of different length !!!")}
  
  lengthNAX <- sum(is.na(X)) # Number of NA values in X
  if(lengthNAX > 0){warning(paste("Vector of true values contains ", lengthNAX, " NA !!! NA excluded", sep = ""))}
  
  lengthNAY <- sum(is.na(Y)) # Number of NA values in Y
  if(lengthNAY > 0){warning(paste("Vector of imputed values contains ", lengthNAY, " NA !!! NA excluded", sep = ""))}
  
  maxX <- max(X, na.rm = T)
  minX <- min(X, na.rm = T)
  distance <- abs(Y-X)/(maxX-minX)
  simi <- 1/(1+distance)
  
  # Check if X is constant and which values of Y are equal to X
  # By definition, replace corresponding 'simi' by 1 (perfect prediction)
  if((maxX-minX)==0){
    warning(print("X is constant"))
    simi[which(Y==mean(X, na.rm = T)&!is.na(X))] <- 1
  }
  
  out <- sum(simi, na.rm = T)/(length(X)-max(lengthNAX, lengthNAY))
  return(out)
  
}