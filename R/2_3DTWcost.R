#' @title DTW-based methods for univariate signals
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description Finds the optimal alignment between two univariate time series based on DTW methods.
#' @param X query vector
#' @param Y response vector
#' @param method "DTW", "DDTW", "AFBDTW", "DTW-D"
#' @param ... additional arguments from functions dtw or dist_afbdtw
#' @import rlist stats dtw
#' @examples
#' data(dataDTWBI)
#' X <- dataDTWBI[, 1] ; Y <- dataDTWBI[, 2]
#' 
#' # Plot query and reference
#' plot(X, type = "l", ylim = c(-5,3))
#' lines(1:length(X), Y, col = "red")
#' 
#' #= Align signals using DTW
#' align_dtw <- minCost(X, Y, method = "DTW")
#' #= Align signals using DDTW
#' align_ddtw <- minCost(X, Y, method = "DDTW")
#' #= Align signals using AFBDTW
#' align_afbdtw <- minCost(X, Y, method = "AFBDTW")
#' #= Align signals using DTW-D
#' align_dtwd <- minCost(X, Y, method = "DTW-D")
#' 
#' #= Plots
#' library(dtw)
#' dtwPlotTwoWay(d = align_dtw, xts <- X, yts = Y, main = "DTW")
#' dtwPlotTwoWay(d = align_ddtw, xts <- X, yts = Y, main = "DDTW")
#' dtwPlotTwoWay(d = align_afbdtw, xts <- X, yts = Y, main = "AFBDTW")
#' dtwPlotTwoWay(d = align_dtwd, xts <- X, yts = Y, main = "DTW-D")
#' 
#' #= Compare cost of each method
#' comparative_cost <- matrix(c(align_dtw$normalizedDistance,
#' align_ddtw$normalizedDistance,
#' align_afbdtw$normalizedDistance,
#' align_dtwd$normalizedDistance), ncol = 4)
#' colnames(comparative_cost) <- c("DTW", "DDTW", "AFBDTW", "DTW-D")
#' comparative_cost

# Cost minimization function
minCost <- function(X, Y, method, ...){
  
  if(length(which(is.na(X)))!=0||length(which(is.na(Y)))!=0){
    stop(print("Input data contain NA"))
  }
  
  method_options <- c("DTW", "DDTW", "AFBDTW", "DTW-D",
                      "dtw", "ddtw", "afbdtw", "dtw-d")
  if(method %in% method_options){
    
    # DTW
  if(method=="DTW"||method=="dtw"){
    out <- dtw(X, Y, ...)
    return(out)
    }
    # DDTW
  if(method=="DDTW"||method=="ddtw"){
    Xtemp <- local.derivative.ddtw(X)
    Ytemp <- local.derivative.ddtw(Y)
    out <- dtw(Xtemp, Ytemp, ...)
    return(out)
    }
    # AFBDTW
  if(method=="AFBDTW"||method=="afbdtw"){
    distTemp <- dist_afbdtw(X, Y, ...)
    out <- dtw(distTemp$distAFBDTW)
    return(out)
    }
    # DTW-D
  if(method=="DTW-D"||method=="dtw-d"){
    out <- dtw(X, Y, keep.internals=TRUE)
    localCostMatrix <- out$localCostMatrix/.ED(X, Y)
    out <- dtw(localCostMatrix)
    return(out)
    }
  }
  
  else(stop("Invalid methodology"))
  
}