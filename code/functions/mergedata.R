#' mergedata : merge expression data using 1:1 orthologs
#' mergedata(a, b, o) --> a_o, b_o, ab_o
mergedata <- function(a, b, o) {
  o_a <- o[, c("one2one", "a")]
  o_b <- o[, c("one2one", "b")]
  merge_x_o <- function(x, o) {
    x_o <- merge(x, o, by.x = 0, by.y = 2)
    rownames(x_o) <- x_o$one2one
    x_o <- x_o[, !sapply(x_o, is.character)]
    x_o <- x_o[order(rownames(x_o)), ]
    return(x_o)
  }
  
  a_o <- merge_x_o(a, o_a)
  b_o <- merge_x_o(b, o_b)
  ab_o <- merge(a_o, b_o, by = 0)
  rownames(ab_o) <- ab_o$Row.names
  ab_o <- ab_o[order(rownames(ab_o)), ]
  ab_o <- ab_o[, !sapply(ab_o, is.character)]
  
  a_o <- a_o[rownames(a_o) %in% rownames(ab_o), ]
  b_o <- b_o[rownames(b_o) %in% rownames(ab_o), ]
  
  res <- list(
    a_o = a_o, 
    b_o = b_o, 
    ab_o = ab_o
  )
  return(res)
}
