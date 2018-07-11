#' @title DTWBI algorithm for univariate signals
#' @author Camille Dezecache, Hong T. T. Phan, Emilie Poisson-Caillault
#' @description Imputes values of a gap of position t_gap and size T in a univariate signal based on DTW algorithm.
#' For more details on the method, see Phan et al. (2017) DOI: <10.1016/j.patrec.2017.08.019>.
#' Default arguments of dtw() function are used but can be manually explicited and modified.
#' @param data input vector containing a large and continuous gap (eventually derived from local.derivative.ddtw() function)
#' @param t_gap location of the begining of the gap (eventually extracted from gapCreation function)
#' @param T_gap gap size (eventually extracted from gapCreation function)
#' @param DTW_method DTW method used for imputation ("DTW", "DDTW", "AFBDTW"). By default "DTW".
#' @param threshold_cos threshold used to define similar sequences to the query. By default, threshold_cos=0.9995 if sequence is longer than 10'000, and threshold_cos=0.995 if shorter.
#' @param thresh_cos_stop Define the lowest cosine threshold acceptable to find a similar window to the query. By default, thresh_cos_stop=0.8.
#' @param step_threshold step used within the loop determining the threshold. By default, step_threshold=50 if sequence is longer than 10'000, step_threshold=10 if sequence length is between 1'000 and 10'000. Else, step_threshold=2.
#' @param ... additional arguments from the dtw() function
#' @return DTWBI_univariate returns a list containing the following elements:
#' \itemize{
#'  \item{output_vector: }{output vector containing complete data including the imputation proposal}
#'  \item{input_vector: }{original vector used as input}
#'  \item{query: }{the query i.e. the adjacent sequence to the gap}
#'  \item{pos_query: }{index of the begining and end of the query}
#'  \item{sim_window: }{vector containing the values of the most similar sequence to the query}
#'  \item{pos_sim_window: }{index of the begining and end of the similar window}
#'  \item{imputation_window: }{vector containing imputed values}
#'  \item{pos_imp_window: }{index of the begining and end of the imputation window}
#' }
#' @import dtw
#' @importFrom lsa cosine
#' 
#' @examples
#' data(dataDTWBI)
#' X <- dataDTWBI[, 1]
#' 
#' rate <- 0.1
#' output <- gapCreation(X, rate)
#' data <- output$output_vector
#' gap_begin <- output$begin_gap
#' gap_size <- output$gap_size

#' imputed_data <- DTWBI_univariate(data, t_gap=gap_begin, T_gap=gap_size)
#' plot(imputed_data$input_vector, type = "l", lwd = 2) # Uncomplete signal
#' lines(imputed_data$output_vector, col = "red") # Imputed signal
#' lines(y = imputed_data$query,
#'       x = imputed_data$pos_query[1]:imputed_data$pos_query[2],
#'       col = "green", lwd = 4) # Query
#' lines(y = imputed_data$sim_window,
#'       x = imputed_data$pos_sim_window[1]:imputed_data$pos_sim_window[2],
#'       col = "orange", lwd = 4) # Similar sequence to the query
#' lines(y = imputed_data$imputation_window,
#'       x = imputed_data$pos_imp_window[1]:imputed_data$pos_imp_window[2],
#'       col = "blue", lwd = 4) # Imputing proposal


DTWBI_univariate <- function(data, t_gap, T_gap, DTW_method = "DTW", threshold_cos = NULL, step_threshold = NULL, thresh_cos_stop = 0.8, ...){

  method_options <- c("DTW", "DDTW", "AFBDTW",
                      "dtw", "ddtw", "afbdtw")
  if(DTW_method %in% method_options){
    print(DTW_method)} else {stop("Invalid DTW method")}
  
  IdGap <- t_gap:(t_gap+T_gap-1)
  if(sum(which(is.na(data[-IdGap])))>0){
    stop("Dataset contains remaining NA outside main gap")
  }
  
  N <- length(data)
  
  if(T_gap>=0.25*N){stop("Gap is to large to compute an appropriate imputation proposal")}
  
  #=== Default parameters definition
  # Default threshold_cos definition
  if(is.null(threshold_cos)){
    if (N>10000) {threshold_cos=0.9995
    } else {threshold_cos=0.995}
  }
  
  # Default step_threshold definition
  if(is.null(step_threshold)){
    if (N>10000) {step_threshold=50
    } else if (N>1000) {step_threshold=10
    } else step_threshold=2
  }
  
  #=== Apply DTW on data ===#
  
  if(t_gap < N/2){ # if t_gap is before the middle of the vector, search right
    
    query_a <- c()
    data_a <- data
    
    # STEP 1: construct a query (temporal window) after the missing data. Q = Dx[t_gap-T_gap;t_gap-1]
    pos_start <- t_gap+T_gap # position of the end of the gap (first value after the gap)
    ind <- pos_start:(pos_start+T_gap-1)
    # Modify query following DTW method used
    query_a <- data[ind]
    if(DTW_method == "DDTW"||DTW_method == "ddtw"){
      query_a <- local.derivative.ddtw(query_a)
      data_a <- local.derivative.ddtw(data_a)
      data_a[t_gap-1] <- 0 # Compensate gap expansion by local.derivative.ddtw() function
      data_a[t_gap+T_gap] <- 0 # Compensate gap expansion by local.derivative.ddtw() function
    }
    
    # STEP 2: find the threshold
    i_start <- pos_start+T_gap
    i_finish <- length(data)-T_gap+1
    
    # selected_qs = indices of selected windows similar to query by features algorithm
    selected_qs <- .DTW_threshold_global_univariate(query = query_a, database = data_a, i_start, i_finish, step_threshold, threshold_cos, thresh_cos_stop, ...)
    
    # STEP 3: find similar windows in the research database		
    if(DTW_method =="AFBDTW"||DTW_method == "afbdtw"){
      id_similar_window <- .Finding_similar_window_univariate_AFBDTW(query = query_a, database = data_a, selected_qs = selected_qs, ...)
    } else {
      id_similar_window <- .Finding_similar_window_univariate(query = query_a, database = data_a, selected_qs = selected_qs, ...)
    }
    
    id_simwin_begin <- id_similar_window[[1]]
    id_simwin_end <- id_simwin_begin+T_gap-1
    similar_query_dtw <- data[id_simwin_begin:id_simwin_end]
    
    # For imputation window
    id_imp_end <- id_simwin_begin-1
    id_imp_begin <- id_imp_end-T_gap+1
    imp_value_dtw <- data[id_imp_begin:id_imp_end]
    
    data_imputed_proposal=data;
    data_imputed_proposal[t_gap:(t_gap+T_gap-1)] <-imp_value_dtw;

    imputation <- list("output_vector" = data_imputed_proposal,
                       "input_vector" = data,
                       "query" = data[ind],
                       "pos_query" = c(pos_start, (pos_start+T_gap-1)),
                       "sim_window" = similar_query_dtw,
                       "pos_sim_window" = c(id_simwin_begin, id_simwin_end),
                       "imputation_window" = imp_value_dtw,
                       "pos_imp_window" = c(id_imp_begin, id_imp_end))
  }
  
  if(t_gap >= N/2){ # if t_gap is after the middle of the vector, search left
    query_b <- c()
    # data_b <- data
    
    # STEP 1: construct a query (temporal window) before the missing data. Q = Dx[t_gap-T_gap;t_gap-1]
    pos_start <- t_gap-1
    ind <- (pos_start-T_gap+1):(pos_start)
    # STEP 2: build a search database before the gap
    Researchbase_b <- data[1:pos_start]
    # Modify query following DTW method used
    query_b <- data[ind]
    if(DTW_method == "DDTW"||DTW_method == "ddtw"){
      query_b <- local.derivative.ddtw(query_b)
      # data_b <- local.derivative.ddtw(data_b)
      # data_b[t_gap-1] <- 0 # Compensate gap expansion by local.derivative.ddtw() function
      # data_b[t_gap+T_gap] <- 0 # Compensate gap expansion by local.derivative.ddtw() function
      Researchbase_b <- local.derivative.ddtw(Researchbase_b)
    }
    
    # STEP 3: find the threshold
    i_start <- 1
    i_finish <- length(Researchbase_b)-2*T_gap # Query excluded from the threshold definition
    # i_finish ends T_gap before the query, so that last reference sequence ends -1 before query
    
    selected_qs <- .DTW_threshold_global_univariate(query = query_b, database = Researchbase_b, i_start=i_start, i_finish=i_finish, step_threshold, threshold_cos, thresh_cos_stop)
    
    # STEP 4: find similar windows in the research database
    if(DTW_method =="AFBDTW"||DTW_method == "afbdtw"){
      id_similar_window <- .Finding_similar_window_univariate_AFBDTW(query = query_b, database = Researchbase_b, selected_qs = selected_qs)
    } else {
      id_similar_window <- .Finding_similar_window_univariate(query = query_b, database = Researchbase_b, selected_qs = selected_qs)
    }
    
    id_simwin_begin <- id_similar_window[[1]]
    id_simwin_end <- id_simwin_begin + T_gap - 1
    similar_query_dtw <- Researchbase_b[id_simwin_begin:id_simwin_end]
    
    # For imputation window
    id_imp_begin <- id_simwin_end + 1
    id_imp_end <- id_imp_begin + T_gap - 1
    imp_value_dtw <- data[id_imp_begin:id_imp_end]
    
    data_imputed_proposal=data;
    data_imputed_proposal[t_gap:(t_gap+T_gap-1)] <-imp_value_dtw;
    
    imputation <- list("output_vector" = data_imputed_proposal,
                       "input_vector" = data,
                       "query" = data[ind],
                       "pos_query" = c((pos_start-T_gap+1), pos_start),
                       "sim_window" = similar_query_dtw,
                       "pos_sim_window" = c(id_simwin_begin, id_simwin_end),
                       "imputation_window" = imp_value_dtw,
                       "pos_imp_window" = c(id_imp_begin, id_imp_end))
  }
  return(imputation)
}