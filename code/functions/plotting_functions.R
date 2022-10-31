### comparABle
### Plotting functions:
#' plot the pca
#' #' pca_sp: plots a pca of the different dataset to see how they group
#' pca_sp(ab_o) --> PCAobject:list(x,y, main responsible genes ...)
plot_pca_ab <- function(pca,ab_o){
  plot(
    pca$x,
    col=c(rep("purple",ncol(ab_o$a_o)),
          rep("green",ncol(ab_o$b_o))
    ),
    pch = 19,
    main = "PCA"
  )
  text(
    pca$x,
    labels=factor(rownames(pca$x)),
    pos=4,
    cex = 0.7
  )
}
