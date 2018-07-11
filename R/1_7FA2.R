#' @title FA2
#' @author Camille Dezecache, Hong T. T. Phan, Emilie Poisson-Caillault
#' @description Estimates the FA2 of two univariate signals Y (imputed values) and X (true values).
#' @details
#' This function returns the value of FA2 of two vectors corresponding to univariate signals X (true values) and Y (imputed values).
#' This FA2 corresponds to the percentage of pairs of values (\eqn{x_{i}, y_{i}}) satisfying the condition \eqn{0,5 <= (Y_{i}/X_{i}) <= 2}.
#' The closer FA2 is to 1, the more accurate is the imputation model.
#' Both vectors Y and X must be of equal length, on the contrary an error will be displayed.
#' In both input vectors, eventual NA will be exluded with a warning diplayed.
#' @param Y vector of imputed values
#' @param X vector of true values
#' @param verbose if TRUE, print advice about the quality of the model
#' @examples
#' data(dataDTWBI)
#' X <- dataDTWBI[, 1] ; Y <- dataDTWBI[, 2]
#' compute.fa2(Y,X)
#' compute.fa2(Y,X, verbose = TRUE)
#' 
#' # By definition, if pairs of true and imputed values are zero,
#' # FA2 corresponding to this pair of values equals 1.
#' X[1] <- 0
#' Y[1] <- 0
#' compute.fa2(Y,X)

compute.fa2<-function(Y, X,verbose=F){
  
  Xtemp <- X
  Ytemp <- Y
  
  if(length(Y)!=length(X)){stop("Input vectors are of different length !!!")}
  
  lengthNAX <- sum(is.na(X)) # Number of NA values
  if(lengthNAX > 0){warning(paste("Vector of true values contains ", lengthNAX, " NA !!! NA excluded", sep = ""))}
  lengthNAY <- sum(is.na(Y)) # Number of NA values
  if(lengthNAY > 0){warning(paste("Vector of imputed values contains ", lengthNAY, " NA !!! NA excluded", sep = ""))}
  
  Xtemp <- X[which(!is.na(X)&!is.na(Y))] # Removing NA
  Ytemp <- Y[which(!is.na(X)&!is.na(Y))] # Removing NA
  
  id0XY <- which(X==0&Y==0) # Identify pairs of XY values equal to 0
  Xtemp[id0XY] <- 1
  Ytemp[id0XY] <- 1
  
  ratio <- (Ytemp/Xtemp)
  fraction <- ratio[ratio>=0.5 & ratio<=2]
  FA2 <- length(fraction)/length(ratio)
  if (verbose){
    if(FA2>0.8) {print("good model")
    }else{print("important number of different points")}
  }
  out<-FA2
  return(out)
}