#' @title Finding similar windows to a query
#' @author Camille Dezecache, Hong T. T. Phan, Emilie Poisson-Caillault
#' @description This function finds similar windows to a query consisting of a univariate vector.
#' @import dtw
#' @importFrom lsa cosine
#' @keywords internal

.Finding_similar_window_univariate <- function(query, database, selected_qs, ...){
  
  # Initialization
  T_gap <- length(query)
  Cost_dist <- c()
  pos_i <- selected_qs
  
  for (i in pos_i){
    k <- i+T_gap-1
    ref <- database[i:k]
    align <- dtw(query, ref, keep.internals=TRUE)
    cost <- align$normalizedDistance
    Cost_dist <- c(Cost_dist, cost)
  }
  
  min_cost <- min(Cost_dist)
  id_similar <- pos_i[which(Cost_dist==min_cost)]
  
  return(list(id_similar))
}