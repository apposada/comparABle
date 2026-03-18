# essentially translate_ids, return.missing turned to FALSE by default.
lookup <- 
  function(x, d, value.col = 2){
    
    d <- d[,c(1,value.col)]
    colnames(d) <- c("key","value")
    d <- d[d$value != "",]
    
    b <- d[match(x,d$key),value.col]
    
    b <- b[complete.cases(b)]
    
    return(b)
    
  }
