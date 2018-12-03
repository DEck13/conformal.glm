

findindex <- function(mat, wn, k){
	A <- seq(from = 0, to = 1 - wn, by = wn)
  n.out <- nrow(mat)
  out <- rep(0, n.out)
  for(j in 1:n.out){
    foo <- 0
    for(i in 1:k){
      if(i < k){
        foo <- foo + 
          (max(which(A < mat[j, i])) - 1) * (1/wn)^(k-i)
      }
      if(i == k){
        foo <- foo + 
          max(which(A < mat[j, i]))
      }
    }
    out[j] <- foo
  }
  out
}



