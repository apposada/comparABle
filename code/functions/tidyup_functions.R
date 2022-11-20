### comparABle
### Tidy-up functions
#' tidyup: remove 0-expressed genes, lowly variable genes etc.

tidyup <- function(x, highlyvariable=F) {
  
  zero_filt <- which(rowSums(x) >= 1 )
  x_zero <- x[zero_filt,]
  cv_x <- apply(x_zero,1,function(z){sd(z)/mean(z)})
  
  if(highlyvariable == T){
    highvar_filt <- which(cv_x > quantile(cv_x, 0.25))
    x_filt <- x_zero[highvar_filt,]
  } else { x_filt <- x_zero}
  
  return(x_filt)
}

#' 
#' rep2means: transform a genexp dataset from rep-based to sample-based
#' basically averages replicates and puts them as one column
#' mirror function of rowMeans_byRepl().
#' colnames must be named 'sample_replicate', and the list 
#' of samples must be named 'sample', in order for the grep
#' step to work.
rowMeans_by_repl <- function(a,genexp){
  numsple <- length(grep(a,colnames(genexp)))
  if(numsple > 1 ){
    c <- rowMeans(genexp[,grep(a, colnames(genexp))])
  } else c <- genexp[,grep(a, colnames(genexp))]
  return(c)
}
rep2means <- function(samples,genexp){
  genexp_mean <- sapply(samples, rowMeans_by_repl,genexp)
  return(genexp_mean)
}

#' 
#' qnorm: prepare a quantile normalisation of input data
#' qnorm(a) --> a_q
#' 
qnorm <- function(x){
  require(preprocessCore)
  colnames_x <- colnames(x)
  rownames_x <- rownames(x)
  x <- as.matrix(x)
  x_q <- normalize.quantiles(x) # figure out why we get zeroes and neg numbers here
  colnames(x_q) <- colnames_x
  rownames(x_q) <- rownames_x
  return(x_q)
}
#' 
#' pair2id --> if no common_name to orthologs pairs,
#' transform a 2-col df of gene_a - gene_b into gene_a - gene_b - pair 
pair2id <- function(x){
  if (ncol(x) == 2){
    numpairs <- nrow(x)
    nDigits <- nchar( trunc( abs(numpairs) ) )
    x$one2one <- paste0(
      "pair_",
      formatC(
        1:numpairs,
        width = nDigits,
        format = "d",
        flag = "0"
      )
    )
    x <- x[,c(3,1,2)]
    colnames(x) <- c("one2one","a","b")
  } else {
    warning("Columns in 1:1 table of orthologs != 2. Did nothing. Please double-check this is correct.")
  }
  return(x)
}
#' 
#' gfam_long --> transform a csv/tsv of OF/OMA/... into long gene-gfam
#' gfam_long(f) --> f:data.frame(gene,gfam)
#' 