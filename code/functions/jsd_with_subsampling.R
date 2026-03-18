#' js_with_subsampling: applies bootstrap using the j-s function
#' 
jsd_with_subsampling <- function(
  a_o,b_o, n = 1000, p = 0.75, bootstrap = FALSE
) {
  
  ensemble_js = vector(mode = "list", length = n)
  
  genes = rownames(a_o)[rownames(a_o) %in% rownames(b_o)]
  
  for (i in 1:n) { 
    if (bootstrap == TRUE) p = 1
    
    ids = sample(genes, p*length(genes), replace = bootstrap)
    
    a_ <- a_o[ids,]
    b_ <- b_o[ids,]
    
    js <- jsd(a_, b_)
    
    ensemble_js[[i]] = js
  }
  
  # calculate the mean matrix
  mean_js <- Reduce("+", ensemble_js) / n
  
  # sd_js <- Reduce("+", lapply(ensemble_js, function(x) (x - apply(ensemble_js, 2, mean))^2)) / n-1
  # sd_js <- sqrt(sd_js)
  
  rownames(mean_js) = colnames(a_o)
  colnames(mean_js) = colnames(b_o)
  
  # rownames(sd_js) = colnames(a_o)
  # colnames(sd_js) = colnames(b_o)
  
  res <- list(
    mean = mean_js#,
    # sd = sd_js
  )
  
  return(res)
  
}
