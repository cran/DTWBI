#' @title Finding similar windows to a query
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description This function finds similar windows to a query consisting of a univariate vector.
#' @import dtw
#' @importFrom lsa cosine
#' @keywords internal

.Finding_similar_window_univariate <- function(query_a, data_a, i_start, i_finish, step_sim_win, threshold, threshold_cos, ...){
  
  T <- length(query_a)
  Dist <- c()
  pos_i <- c()
  
  id_found_start <- c()
  Cosine_threshold <- c()
  CC <- c()
  dist_found <- c()
  i <- i_start
  query_gf <- .globalfeatures(query_a)
  
  while(i < i_finish){
    k <- i+T-1
    ref <- data_a[i:k]
    lag <- round(0.1*length(ref), 0)
    cc <- stats::ccf(query_a, ref, lag.max=lag, plot=F)$acf
    cc <- max(abs(cc))
    CC <- c(CC,cc)
    
    ref_gf <- .globalfeatures(ref)
    
    if((T<=1000)||all(query_a==0)||all(ref==0)){
      align <- dtw(query_a, ref, ...)
      cost <- align$normalizedDistance
      Dist <- c(Dist, cost)
      if (cost<=threshold){
        id_found_start <- c(id_found_start, i)
        dist_found <- c(dist_found, cost)
      }
    } else {
      cos_threshold <- cosine(as.numeric(query_gf), as.numeric(ref_gf))
      Cosine_threshold <- c(Cosine_threshold, cos_threshold[1])
      pos_i <- c(pos_i, i)
      if(abs(cos_threshold[1])>=threshold_cos){
        align <- dtw(query_a, ref, ...)
        cost <- align$normalizedDistance
        Dist <- c(Dist, cost)
        if(cost<=threshold){
          id_found_start <- c(id_found_start, i)
          dist_found <- c(dist_found, cost)
        }
      }
    }
    i <- i+step_sim_win
  }
  
  Cc <- CC
  
  if(length(id_found_start)<=30){
    return(list(id_found_start, dist_found, Cc))
  }else{
    id_similar <- sample(id_found_start, 30, replace=FALSE)
    p1 <- c()
    for (i in 1:30){
      t1 <- which(id_found_start==id_similar[i])
      p1 <- c(p1, dist_found[t1])
    }
    return (list(id_similar, p1, Cc))
  }
}