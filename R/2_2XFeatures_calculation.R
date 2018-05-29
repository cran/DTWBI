# Local feature
.local_feature <- function(q){
  n <- length(q)
  out <- matrix(nrow = n, ncol = 2)
  out[1, ] <- q[1]
  for (i in 2:(n-1))
  {
    t <- c(q[i]-q[(i-1)], q[i]-q[(i+1)])
    out[i, ] <- t
  }
  out[n, ] <- q[n]
  return(out)
}

# Global feature
.global_feature <- function(q){
  n <- length(q)
  out <- matrix(nrow = n, ncol = 2)
  out[1, ] <- q[1]
  for (i in 2:(n-1))
  {
    t <- c(q[i]-mean(q[1:(i-1)], na.rm = T),q[i]-mean(q[(i+1):n], na.rm = T))
    out[i, ] <- t
  }
  out[n, ] <- q[n]
  return(out)
}

# Distance matrix
.manhattan_dist <- function(p, q){
  return(sum(abs(p-q), na.rm = T))
}

.dist_matrix <- function(q, r){
  m <- nrow(q)
  n <- nrow(r)
  dist <- matrix(10000, m, n)
  
  for(i in 2:(m-1)){
    for (j in 2:(n-1)){
      t <- .manhattan_dist(q[[i]], r[[j]])
      dist[i,j] <- t
      
      # dist[i, j] <- abs(qfeature[i, 1]-rfeature[i, 1]) + abs(qfeature[i, 2]-rfeature[i, 2])
    }
  }
  return(dist)
}