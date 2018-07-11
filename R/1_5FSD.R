#' @title Fraction of Standard Deviation (FSD)
#' @author Camille Dezecache, Hong T. T. Phan, Emilie Poisson-Caillault
#' @description Estimates the Fraction of Standard Deviation (FSD) of two univariate signals Y (imputed values) and X (true values).
#' @details
#' This function returns the value of FSD of two vectors corresponding to univariate signals.
#' Values of FSD closer to zero indicate a better performance method for the imputation task.
#' Both vectors Y and X must be of equal length, on the contrary an error will be displayed.
#' In both input vectors, eventual NA will be exluded with a warning diplayed.
#' @param Y vector of imputed values
#' @param X vector of true values
#' @param verbose if TRUE, print advice about the quality of the model
#' @importFrom stats sd
#' @examples
#' data(dataDTWBI)
#' X <- dataDTWBI[, 1] ; Y <- dataDTWBI[, 2]
#' compute.fsd(Y,X)
#' compute.fsd(Y,X, verbose = TRUE)
#' 
#' # By definition, if true and imputed values are equal and constant,
#' # FSD = 0.
#' X <- rep(runif(1), 10)
#' Y <- X
#' compute.fsd(Y,X)
#' 
#' # However, if true and imputed values are constant but different,
#' # FSD is not calculable. An error is displayed.
#' \dontrun{
#' X <- rep(runif(1), 10);Y <- rep(runif(1), 10)
#' compute.fsd(Y,X)}

compute.fsd <- function(Y, X, verbose=F){
  
  if(length(Y)!=length(X)){stop("Input vectors are of different length !!!")}
  
  lengthNAX <- sum(is.na(X)) # Number of NA values
  if(lengthNAX > 0){warning(paste("Vector of true values contains ", lengthNAX, " NA !!! NA excluded", sep = ""))}
  lengthNAY <- sum(is.na(Y)) # Number of NA values
  if(lengthNAY > 0){warning(paste("Vector of imputed values contains ", lengthNAY, " NA !!! NA excluded", sep = ""))}
  
  sd1=sd(Y, na.rm = T)
  sd2=sd(X, na.rm = T)
  
  if(sd1==sd2){
    if(mean(X, na.rm = T)==mean(Y, na.rm = T)){
      warning("Vectors of true and imputed values are constant and equal !!! By definition FSD=0")
      FS <- 0
      if(verbose){print("acceptable model")}
      }
    else{
      stop("Impossible to compute FSD: vectors of true and imputed values are constant but different !!!")
    }
    out <- FS
    return(out)
  }
  else{
  FS <- 2*abs((sd1-sd2)/(sd1+sd2))
  if(verbose){
    if(abs(FS)<0.5) {print("acceptable model");
    }else{print("non acceptable FS");}
  }
  out <- FS
  return(out)
  }
}