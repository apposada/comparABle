### comparABle
### Basic calculation functions
#' jsd: perform jensen-shannon distance; MUST FIND BETTER WAY
jsd <- function(a,b){
  
  js <- function(p,q){ #jensen-shannon from https://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r
    n <- 0.5 * (p + q)
    JS <- sqrt(0.5 * (sum(p * log((p / n)+0.1)   + sum(q * (log((q / n)+0.1)))))) #sqrt(0.5 * (sum(p * log(p / n)) + sum(q * log(q / n)))) # with the original one I get NaNs in some, because I get log(0), log(-1), etc.
    return(JS)
  }
  
  mat <- matrix(NA,nrow=ncol(a),ncol=ncol(b))
  for (i in 1:ncol(a)){
    for (j in 1:ncol(b)){
      mat[i,j] <- js(a[,i],b[,j])
    }
  }
  return(mat)
}