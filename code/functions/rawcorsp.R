#' rawcorsp : calculate correlation/distance between spp 
#' using several or user-defined metrics (default: pearson, 
#' spearman, jensen-shannon distance)
#' rawcorsp(a_o, b_o) --> list(pe:matrix, sp:matrix, js:matrix)
rawcorsp <- function(a_o, b_o){ # make one for TFs or genes of interest
  pe <- cor(a_o, b_o, m = "pe")
  sp <- cor(a_o, b_o, m = "sp")
  js <- jsd(a_o, b_o)
  
  res <- list(
    pe = pe, 
    sp = sp, 
    js = js
  )
  return(res)
}
