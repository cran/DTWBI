#' @title Global threshold for missing data imputation
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description Finds a threshold for univariate missing data imputation in a univariate vector.
#' @import dtw
#' @importFrom lsa cosine
#' @keywords internal

.DTW_threshold_global_univariate_AFBDTW <- function(query_a, data_a, i_start, i_finish, step_threshold, threshold_cos, ...){

  Dist_global_threshold <- c()
  T <- length(query_a)
  
  i <- i_start
  Cosine_threshold <- c()
  pos_i <- c()
  
  while(i < i_finish){ 
    k <- i + T - 1
    ref <- data_a[i:k]	
    if((T < 600)||all(query_a==0)||all(ref==0)) 
    {
      dist_matrix_AFBDTW <- dist_afbdtw(query_a, ref)
      align <- dtw(dist_matrix_AFBDTW$distAFBDTW, ...)
      cost <- align$normalizedDistance
      Dist_global_threshold <- c(Dist_global_threshold, cost)
    }
    else{
      cos_threshold <- cosine(as.numeric(.globalfeatures(query_a)), as.numeric(.globalfeatures(ref)))
      Cosine_threshold <- c(Cosine_threshold, cos_threshold[1])
      pos_i <- c(pos_i, i)	
      if(abs(cos_threshold[1])>=threshold_cos){
        dist_matrix_AFBDTW <- dist_afbdtw(query_a, ref)
        align <- dtw(dist_matrix_AFBDTW$distAFBDTW, ...)
        cost <- align$normalizedDistance
        Dist_global_threshold <- c(Dist_global_threshold, cost)
      }
    }
    i <- i + step_threshold
  }
  while(length(Dist_global_threshold)==0){
    threshold_cos <- threshold_cos-0.01
    len <- length(pos_i)
    for (j in 1:len){
      if(abs(Cosine_threshold[j])>=threshold_cos){
        i <- pos_i[j]
        k <- i + T -1
        ref <- data_a[i:k]
        dist_matrix_AFBDTW <- dist_afbdtw(query_a, ref)
        align <- dtw(dist_matrix_AFBDTW$distAFBDTW, ...)
        cost <- align$normalizedDistance
        Dist_global_threshold <- c(Dist_global_threshold, cost)
      }
    }
  }
  
  global_threshold <- min(Dist_global_threshold)
  return(list((global_threshold), threshold_cos))
}