#' @title Global threshold for missing data imputation
#' @author Camille Dezecache, Hong T. T. Phan, Emilie Poisson-Caillault
#' @description Finds a threshold for univariate missing data imputation in a univariate vector.
#' @import dtw
#' @importFrom lsa cosine
#' @keywords internal

.DTW_threshold_global_univariate <- function(query, database, i_start, i_finish, step_threshold, threshold_cos, thresh_cos_stop, ...){
  
  # Initialization
  T_gap <- length(query)
  Cosine_threshold <- c()
  pos_i <- c()
  threshold_cos_temp <- threshold_cos
  
  while((length(Cosine_threshold)==0)&&(threshold_cos_temp>thresh_cos_stop)){
    i <- i_start
    while(i<=i_finish){ 
      k <- i+T_gap-1
      ref <- database[i:k]
      
      gf_q <- as.numeric(.globalfeatures(query))
      gf_r <- as.numeric(.globalfeatures(ref))
      
      ind_nan <- NULL
      ind_nan <- c(which(is.nan(gf_q)), which(is.nan(gf_r))) # Remove NaN in global features
      
      if(length(ind_nan)>0){
        gf_q <- gf_q[-ind_nan];
        gf_r <- gf_r[-ind_nan]
      }
      
      cos_threshold <- cosine(gf_q, gf_r)
      
      if(cos_threshold[1]>=threshold_cos_temp){
        pos_i <- c(pos_i, i)
        Cosine_threshold <- c(Cosine_threshold, cos_threshold[1])
      }
      i <- i+step_threshold
    }
    threshold_cos_temp <- threshold_cos_temp-0.01
  }
  
  if(length(Cosine_threshold)==0){stop("No similar window looks appropriate for imputation")}
  
  return(pos_i)
}