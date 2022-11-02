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

#' Functions for heatmaps of species A and B, ideally a color annotation of
#' purple and green for the correlations, the orthology strategy, etc.
#' 
#' Functions for gene age enrichment (probably change the colours used in
#' that analysis from purplegreen to something else to avoid overlap with
#' species notation)
#' 
#' Functions for gene age barplot for each species
#' 
#' Functions for plotting of COGs per species; bot barplots, heatmaps, pies