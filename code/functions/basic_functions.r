### comparABle
### Basic calculation functions
#' jsd: perform jensen-shannon distance; MUST FIND BETTER WAY
jsd <- function(a,b){
  
  a <- apply(a+1,2,function(x){x/sum(x)})
  b <- apply(b+1,2,function(x){x/sum(x)})
  
  js <- function(p,q){ #jensen-shannon from https://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r
    n <- 0.5 * (p + q)
    JS <- sqrt(0.5 * (sum(p * log((p / n)) + sum(q * (log((q / n)))))))
    return(JS)
  }
  
  mat <- matrix(NA,nrow=ncol(a),ncol=ncol(b))
  for (i in 1:ncol(a)){
    for (j in 1:ncol(b)){
      mat[i,j] <- js(a[,i],b[,j])
    }
  }
  
  rownames(mat) <- colnames(a)
  colnames(mat) <- colnames(b)
  
  return(mat)
}
