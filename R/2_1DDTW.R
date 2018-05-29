#' @title Local derivative estimate to compute DDTW
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description This function estimates the local derivative of a vector.
#' It can be used as an input in dtw() function (package dtw) to align two univariate signals.
#' @param X input vector from which local derivative has to be calculated
#' @examples
#' data(dataDTWBI)
#' X <- dataDTWBI[, 1]
#' local.derivative.ddtw(X)
#' 
#' # Plot
#' plot(X, type = "b", ylim = c(-1, 1))
#' lines(local.derivative.ddtw(X), col = "red")

local.derivative.ddtw <- function(X){
  out <- X
  for (i in 2:(length(X)-1))
  {
    tmp <- ((X[i]-X[i-1])+(X[i+1]-X[i-1])/2)/2
    out[i] <-  tmp
  }
  return(out)
}