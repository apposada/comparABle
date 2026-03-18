#' @param df an orthogroup x species data frame containing at least your two species of interest, with gene ids separated as comma-separated
#' @param sp_a name/ID of the species 'a'
#' @param sp_b name/ID of the species 'b'
#' @returns a gene x orthogroup data frame of both species, long-format, non-comma separated.
fam_longformat <- function(df, sp_a, sp_b){
  
  # a
  col_a = which(colnames(df) == sp_a)
  f_a <- 
    df[df[,col_a] != "",c(1,col_a)]
  f_a = 
    stack(lapply(
      split(f_a[[sp_a]],f_a[,1]), # a temporary list resulting from splitting all values in "sp_a" column by OG id
      function(x){strsplit(x,",")[[1]]} # within each element of that temporary list, we split all comma-separated elements
    ))#[,c(2,1)] # the stack above puts this back in place as a long, "melted" format
  
  # b # do the same for species b
  col_b = which(colnames(df) == sp_b)
  f_b <- 
    df[df[,col_b] != "",c(1,col_b)]
  f_b = 
    stack(lapply(
      split(f_b[[sp_b]],f_b[,1]),
      function(x){strsplit(x,",")[[1]]}
    ))#[,c(2,1)]
  
  # put together, rename tidy
  Y = rbind(f_a,f_b)
  colnames(Y) <- c("id","gfam")
  Y$gfam <- as.character(Y$gfam)
  
  return(Y)
}
