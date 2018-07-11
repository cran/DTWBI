#' @title Gap creation
#' @author Camille Dezecache, Hong T. T. Phan, Emilie Poisson-Caillault
#' @description This function creates a large continuous gap within a univariate signal.
#' Gap size is defined as a percentage of input vector length.
#' By default, the created gap starts at a random location.
#' @param X input vector
#' @param rate size of desired gap, as a percentage of input vector size
#' @param begin location of the begining of the gap (random by default)
#' @return gapCreation returns a list containing the following elements:
#' \itemize{
#'  \item{output_vector: }{output vector containing the created gap}
#'  \item{input_vector: }{original vector used as input}
#'  \item{begin_gap: }{index of the begining of the gap}
#'  \item{rate: }{size of the created gap in percentage of the input vector length}
#'  \item{gap_size: }{length of the created gap}
#' }
#' @examples
#' data(dataDTWBI)
#' X <- dataDTWBI[, 1]
#' rate <- 0.1
#' output <- gapCreation(X, rate)
#' plot(output$input_vector, type = "l", col = "red", lwd = 2)
#' lines(output$output_vector, lty = "dashed", lwd = 2)

gapCreation <- function(X, rate, begin=NULL){
  
  Xgap <- X
  
  # No missing data (pass-through)
  if(rate == 0){
    warning("rate = 0, no gap created")
    return(X)
  }
  
  if(rate == 1){
    stop("rate must be < 1")
  }
  
  gap_size <- round(rate*length(X))
  if(is.null(begin)){
    gap_id <- sample(1:(length(X)-gap_size), 1)
  } else {gap_id <- begin}
  Xgap[gap_id:(gap_id+gap_size-1)] <- NA
  
  gap_final <- list("output_vector" = Xgap,
                    "input_vector" = X,
                    "begin_gap" = gap_id,
                    "rate" = rate,
                    "gap_size" = gap_size)
  return(gap_final)
}