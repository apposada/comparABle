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


#' plot_hicor_genes:
plot_hicor_genes <- function(scatter,GO_table,age_table,age_enrichemnt_hm){
  require(ggplot2)
  require(cowplot)
  require(RColorBrewer)
  require(ComplexHeatmap)
  require(readr)
  require(reshape2)

  # create barplot for GO terms
  colormap <- colorRampPalette(viridis::magma(11))
  GO_table$classicFisher <- parse_number(GO_table$classicFisher)
  GO_plot <-
    ggplot(
      GO_table, 
      aes(
        x = reorder(Term, -log10(classicFisher)),
        y = Significant/Expected,
        fill = -log10(classicFisher)
        )
      )+
    geom_bar(stat = "identity") +
    scale_fill_gradientn(colors = colormap(100), name = "-log10(pvalue)") +
    coord_flip() +
    theme_classic() +
    labs(x = "GO term", y = "Significant/Expected ratio")
  
  # create barplot for age data
  age_plot <-
    ggplot(reshape2::melt(age_table[1,]), aes(x = variable, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "Age", y = "Number of genes")+
    scale_fill_brewer(
      type = "div",
      palette = "Spectral",
      direction = 1
    )
  # create heatmap drawing
  hm <- grid.grabExpr(draw(age_enrichemnt_hm))
  
  # combine all plots into a grid
  plot_grid(
    plot_grid(scatter, GO_plot,ncol = 2, align = "v", axis = "tb", rel_widths = c(1,2.5)),
    plot_grid(age_plot, hm, ncol = 2, align = "hv", axis = "tb"),
    ncol = 1
  )
}


#' 
plot_cors <- function(cors){
  h1 <- 
    Heatmap(
      cors$pe,
      name = "Pearson",
      cluster_rows = F, cluster_columns = F,
      show_row_names = TRUE,
      col = sequential_hcl(10,"BluYl", rev = TRUE)
      )
  
  h2 <- 
    Heatmap(
      cors$sp,
      name = "Spearman",
      cluster_rows = F, cluster_columns = F,
      show_row_names = TRUE,
      col = sequential_hcl(10,"YlOrRd", rev = TRUE)
      )
  
  h3 <- 
    Heatmap(
      cors$js,
      name = "JSD", 
      cluster_rows = F, cluster_columns = F, 
      show_row_names = TRUE, 
      col = brewer.pal(10,"RdBu")
      )
  
  h_list <- h1+h2+h3
  
  draw(h_list, auto_adjust = FALSE)
}
