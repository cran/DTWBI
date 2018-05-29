#' @title Fractional Bias (FB)
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description Estimates the Fractional Bias (FB) of two univariate signals Y (imputed values) and X (true values).
#' @details
#' This function returns the value of FB of two vectors corresponding to univariate signals, indicating whether predicted values are underestimated or overestimated compared to true values.
#' A perfect imputation model gets \eqn{FB = 0}.
#' An acceptable imputation model gives \eqn{FB <= 0.3}.
#' Both vectors Y and X must be of equal length, on the contrary an error will be displayed.
#' In both input vectors, eventual NA will be exluded with a warning diplayed.
#' @param Y vector of imputed values
#' @param X vector of true values
#' @param verbose if TRUE, print advice about the quality of the model
#' @examples
#' data(dataDTWBI)
#' X <- dataDTWBI[, 1] ; Y <- dataDTWBI[, 2]
#' compute.fb(Y,X)
#' compute.fb(Y,X, verbose = TRUE)
#' 
#' # If mean(X)=mean(Y)=0, it is impossible to estimate FB,
#' # unless both true and imputed values vectors are constant.
#' # By definition, in this case, FB = 0.
#' X <- rep(0, 10) ; Y <- rep(0, 10)
#' compute.fb(Y,X)
#' 
#' # If true and imputed values are not zero and are opposed, FB = Inf.
#' X <- rep(runif(1), 10)
#' Y <- -X
#' compute.fb(Y,X)

compute.fb<-function(Y, X,verbose=F){
  
  if(length(Y)!=length(X)){stop("Input vectors are of different length !!!")}
  
  lengthNAX <- sum(is.na(X)) # Number of NA values
  if(lengthNAX > 0){warning(paste("Vector of true values contains ", lengthNAX, " NA !!! NA excluded", sep = ""))}
  lengthNAY <- sum(is.na(Y)) # Number of NA values
  if(lengthNAY > 0){warning(paste("Vector of imputed values contains ", lengthNAY, " NA !!! NA excluded", sep = ""))}
  
  m1=mean(Y, na.rm= T)
  m2=mean(X, na.rm = T)
  
  if(m1!=0&&m1==-m2){warning("X=-Y => FB=Inf")}
  
  if(m1==0&&m2==0){
    if((max(Y, na.rm = T)-min(Y, na.rm = T))==(max(X, na.rm = T)-min(X, na.rm = T))){
      warning("Vectors of input and imputed values are equal and constant. By definition FB=0")
      FB <- 0
      if(verbose){print("acceptable FB")}
      out<-FB
    }
    else{stop(print("Impossible to estimate FB: vectors of input and imputed values have the same mean 0 but are not constant !!!"))}
  }
  else{
    FB <- 2*abs((m1-m2)/(m1+m2))
    if (verbose){
      if(abs(FB)<0.3) {print("acceptable FB")
      }else{print("non acceptable FB")}
    }
    out<- FB
  }
  
  return(out)
}