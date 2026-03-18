#' get_high_cor_genes: function to retrieve highly expressed genes
#' across top-similar pairs of stages. Sub-optimal;
#' workaround until we manage to implement treegapgenes function
get_high_cor_genes <-
  function(
    mat, a_o,b_o, o = NULL, weights_method = "neg_exp",
    topgenes_filt_method = "lm", top_N = 5, stages = NULL, ...
  ) {
    require(ggpointdensity)
    require(ggplot2)
    
    mat <- mat
    ab_o <- cbind(a_o,b_o)
    
    
    if (is.null(stages)) {
    # Find the indices of the top N highest elements
    ind <- head(order(mat, decreasing = FALSE), top_N)
    
    # Convert the indices to rows and columns
    row_ind <- (ind - 1) %% nrow(mat) + 1
    col_ind <- (ind - 1) %/% nrow(mat) + 1
    
    # Print the result
    cat("The top five similar pair of stages are:", mat[ind], "\n")
    } else {
      if (length(stages$a) != length(stages$b)) stop("incorrect number of stages chosen, are you missing something?")
      row_ind <- sapply(stages$a,function(x){grep(x,colnames(a_o))})
      col_ind <- sapply(stages$b,function(x){grep(x,colnames(b_o))})
      ind <- 1:length(stages$a)
    }
    
    print("Subsetting count matrices")
    hcor_ab <- list()
    for (i in 1:length(ind)) {
      
      hco <-
        data.frame(
          a = log1p(a_o[,row_ind[i]]),
          b = log1p(b_o[,col_ind[i]]),
          row.names = rownames(a_o)
        )
      
      hcor_ab[[i]] <- hco
      
      names(hcor_ab)[i] <- 
        paste0(
          "cor_",
          colnames(a_o)[row_ind[i]],
          "__",
          colnames(b_o)[col_ind[i]]
        )
    }
    
    print("Model fitting and top genes")
    hcor_ab_topgenes <- list()
    for (i in 1:length(hcor_ab)){
      
      jsd_value <- round(mat[row_ind[i],col_ind[i]],2)
      
      hco <- hcor_ab[[i]]
      
      name_ <- names(hcor_ab)[i]
      
      # lm
      hco_lm <- lm(b ~ a, data = hco)
      
      #define weights to use
      
      if (weights_method == "neg_exp") {
        # use relative of neg exponential of square fitted values as weights
        wt <- exp(-lm(abs(hco_lm$residuals) ~ hco_lm$fitted.values)$fitted.values)
        wt <- relativise(wt)
      } else {
        # else use the square fit of the relationship between residuals and fits
        wt <- 1 / lm(abs(hco_lm$residuals) ~ hco_lm$fitted.values)$fitted.values^2
      }
      
      #perform weighted least squares regression
      hco_lm_wls <- lm(b ~ a, data = hco, weights = wt)
      
      if (topgenes_filt_method == "lm"){
        
        colnames_of_interest <- c(
          colnames(a_o)[row_ind[i]],
          colnames(b_o)[col_ind[i]]
        )
        
        topgenes <- topgenes_lmfit(
          data = ab_o,
          cols_of_interest = colnames_of_interest
        )
        
        filt <- which(rownames(hco) %in% topgenes)
        
      } else if (topgenes_filt_method == "rank"){
        #rank, will use later
        hco$sum <- hco$a + hco$b
        hco$rank <- rank(-hco$sum,ties.method = "min")
        filt <- which(hco$rank <= quantile(hco$rank,.1))
      } else {
        filt <- which(
          abs(hco_lm_wls$residuals) < 1 &
            leverage > mean(leverage) &
            hco$a > mean(hco$a) &
            hco$b > mean(hco$b) 
        )
      }
      
      # plots
      {
        ab_scatter_ggplot2 <-
          ggplot(data = hco,mapping = aes(x=a,y=b)) +
          labs(title = name_ ) +
          geom_pointdensity()+
          stat_smooth(method = "lm", se = FALSE, color = "red") + 
          theme_minimal()+
          scale_color_viridis()+
          annotate(
            "text", x = Inf, y = -Inf, hjust = 1, vjust = 0, 
            label = paste(
              "Equation: y =", 
              sprintf("%.2f", coef(summary(lm(b ~ a, data = hco, weights = wt)))[2, 1]), "x +", 
              sprintf("%.2f", coef(summary(lm(b ~ a, data = hco, weights = wt)))[1, 1]), "\n",
              "J-S Divergence = ", jsd_value
            )
          )
        
        p <- 
          ggplot(data = hco, aes(x = a, y = b)) +
          geom_point(pch = 20, col = rgb(0.1,0.1,0.1,0.1)) +
          geom_abline(
            slope = hco_lm$coefficients[2], 
            intercept = hco_lm$coefficients[1], col = "blue") +
          geom_abline(
            slope = hco_lm_wls$coefficients[2], 
            intercept = hco_lm_wls$coefficients[1], col = "red", lwd = 1.25) +
          geom_point(data = hco[filt,], pch = 21, col = "black", bg = "lightgreen") +
          labs(title = name_)+
          theme_minimal()
      }
      
      if (is.null(o) == TRUE) {
        top_genes <- rownames(hco)[filt]
      } else{
        top_genes <- o[o$one2one %in% rownames(hco)[filt],]
      }
      
      res_i <- list(
        top_genes = top_genes,
        lm = hco_lm_wls,
        plot = ab_scatter_ggplot2,
        plot_topgenes = p
      )
      
      hcor_ab_topgenes[[i]] <- res_i
      names(hcor_ab_topgenes)[i] <- name_
    }
    
    res <- list(
      hicor_matrices = hcor_ab,
      hicor_topgenes = hcor_ab_topgenes
    )
    
    return(res)
  }
