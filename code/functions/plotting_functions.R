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

plot_hicor_genes <- function(scatter,GO_table,age_table,age_enrichemnt_hm,GO_plot = NULL){
  require(ggplot2)
  require(cowplot)
  require(RColorBrewer)
  require(ComplexHeatmap)
  require(readr)
  require(reshape2)
  
  if (is.null(GO_plot)) {
    # create barplot for GO terms
    colormap <- colorRampPalette(viridis::magma(11))
    if (class(GO_table$classicFisher) == "character") {
      GO_table$classicFisher <- parse_number(GO_table$classicFisher)
    }
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
  } else{
    GO_plot <- GO_plot
  }
  
  
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

require(RColorBrewer)
plot_cors <- function(cors){
  require(RColorBrewer)
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

hm_js <-
  function(m){
    h <- 
      Heatmap(
        m, 
        col = brewer.pal(10,"RdBu"),
        cluster_rows = F,
        cluster_columns = F,
        name = "JSD",
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6)
      )
    h
  }

hm_hypg <-
  function(m){
    h <-
      Heatmap(
        m,
        col = sequential_hcl(10,"YlOrRd", rev = TRUE),
        cluster_rows = F,
        cluster_columns = F,
        name = "-log(pval)",
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6)
      )
    h
  }

# plot_ab_scatter_ggplot2 <- function(){
#   require(ggplot2)
#   ggplot(data = hco,mapping = aes(x=a,y=b)) +
#     labs(title = name_ ) +
#     geom_pointdensity()+
#     stat_smooth(method = "lm", se = FALSE, color = "red") + 
#     theme_minimal()+
#     scale_color_viridis()+
#     annotate(
#       "text", x = Inf, y = -Inf, hjust = 1, vjust = 0, 
#       label = paste(
#         "Equation: y =", 
#         sprintf("%.2f", coef(summary(lm(b ~ a, data = hco, weights = wt)))[2, 1]), "x +", 
#         sprintf("%.2f", coef(summary(lm(b ~ a, data = hco, weights = wt)))[1, 1]), "\n",
#         "J-S Divergence = ", jsd_value
#       )
#     )
# }

# plot_ab_scatter_topgenes <- function(){
#   require(ggplot2)
#   ggplot(data = hco, aes(x = a, y = b)) +
#     geom_point(pch = 20, col = rgb(0.1,0.1,0.1,0.1)) +
#     geom_abline(
#       slope = hco_lm$coefficients[2], 
#       intercept = hco_lm$coefficients[1], col = "blue") +
#     geom_abline(
#       slope = hco_lm_wls$coefficients[2], 
#       intercept = hco_lm_wls$coefficients[1], col = "red", lwd = 1.25) +
#     geom_point(data = hco[filt,], pch = 21, col = "black", bg = "lightgreen") +
#     labs(title = name_)+
#     theme_minimal()
# }

