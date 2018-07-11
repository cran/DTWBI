#' @title Adaptive Feature Based Dynamic Time Warping algorithm
#' @author Camille Dezecache, Hong T. T. Phan, Emilie Poisson-Caillault
#' @description 
#' This function estimates a distance matrix which is used as an input in dtw() function (package dtw) to align two univariate signals following Adaptative Feature Based Dynamic Time Warping algorithm (AFBDTW).
#' @param q query vector
#' @param r reference vector
#' @param w1 weight of local feature VS global feature.
#'  By default, w1 = 0.5, and by definition, w2 = 1 - w1.
#' @return A list containing the following elements:
#' \itemize{
#' \item{query: }{the query vector}
#' \item{response: }{the response vector}
#' \item{query_local: }{local feature of the query}
#' \item{response_local: }{local feature of the response vector}
#' \item{query_global: }{global feature of the query}
#' \item{response_global: }{global feature of the response vector}
#' \item{dist_local: }{distance matrix of the local feature}
#' \item{dist_local: }{distance matrix of the global feature}
#' \item{distAFBDTW: }{AFBDTW distance matrix}
#' }
#' @import rlist stats
#' @examples
#' data(dataDTWBI)
#' X <- dataDTWBI[, 1] ; Y <- dataDTWBI[, 2]
#' AFBDTW_Dist <- dist_afbdtw(X, Y)

dist_afbdtw <- function(q, r, w1=0.5){
  
  w2 <- 1-w1
  
  if(w1<=0){stop("Weights should be positive")}
  
  ql <- .local_feature(q)
  rl <- .local_feature(r)
  qg <- .global_feature(q)
  rg <- .global_feature(r)
  dist_local <- .dist_matrix(ql, rl)
  dist_global <- .dist_matrix(qg, rg)
  dist <- w1*dist_local + w2*dist_global
  
  outputAFBDTW <- list("query" = q,
                       "response" = r,
                       "query_local" = ql,
                       "response_local" = rl,
                       "query_global" = qg,
                       "response_global" = rg,
                       "dist_local" = dist_local,
                       "dist_global" = dist_global,
                       "distAFBDTW" = dist)
}