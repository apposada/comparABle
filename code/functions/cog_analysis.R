### comparABle
### COG analysis

require(ComplexHeatmap)
require(circlize)
require(ggplot2)
require(dplyr)

#' function to calculate cog enrichment for list  of N group of genes x J cogs.
#' 
#' @param x_modules 
#' 
cog_enrichment_analysis <- function(x_modules, x_cog, specific_cogs = FALSE, fisher_pval = 0.05, plotting = FALSE, functional_categories = FALSE) {
    modules <- sort(unique(x_modules$module))
    if (specific_cogs != FALSE) {
        cogs <- sort(unique(specific_cogs)) # must have numerics at the beginning
    } else {
        cogs <- sort(unique(x_cog$cog)) # must have numerics at the beginning
    }

    x_module_cog <- unique(
        merge(
            x_modules,
            x_cog,
            by = 1
        )
    )

    quantmat <- matrix(
        rep(0,times = length(modules)*length(cogs)),
        nrow = length(modules),
        ncol = length(cogs)
    )

    enrichmat <- matrix(
        rep(1,times = length(modules)*length(cogs)),
        nrow = length(modules),
        ncol = length(cogs)
    )

    pvalmat <- matrix(
        rep(1,times = length(modules)*length(cogs)),
        nrow = length(modules),
        ncol = length(cogs)
    )

    for (i in modules) {
        h <- which(modules == i)
        for (j in cogs) {
            n <- which(cogs == j)

            contingency <- matrix(
                c(
                    sum(x_module_cog$module == i & x_module_cog$cog == j),
                    sum(x_module_cog$module == i & x_module_cog$cog != j),
                    sum(x_module_cog$module != i & x_module_cog$cog == j),
                    sum(x_module_cog$module != i & x_module_cog$cog != j)
                    ),
                nrow = 2
            )

            quantmat[h,n] <- contingency[1,1] # of cog j genes in module i

            enrichmat[h,n] <- 100 * # percent
            (
                contingency[1,1]/(contingency[2,1] + contingency[1,1]) -
                contingency[1,2]/(contingency[2,2] + contingency[1,2])
            )

            pvalmat[h,n] <- fisher.test(contingency)$p.value
        }
    }

    quantdf <- as.data.frame(quantmat)
    colnames(quantdf) <- cogs; rownames(quantdf) <-  modules

    enrichdf <- as.data.frame(enrichmat)
    colnames(enrichdf) <- cogs; rownames(enrichdf) <-  modules

    pvaldf <- as.data.frame(pvalmat)
    colnames(pvaldf) <- cogs ; rownames(pvaldf) <- modules
    pvaldf[pvaldf > fisher_pval ] <- 1 # clipped
    
    if(functional_categories == TRUE){
        cogs_functionalcategories <- cogs_functionalcategories[
            cogs_functionalcategories$cog %in% cogs,
        ]
        colnames(quantdf) <- cogs_functionalcategories$functional_category
        colnames(enrichdf) <- cogs_functionalcategories$functional_category
        colnames(pvaldf) <- cogs_functionalcategories$functional_category
    }

    col_cogs <- colorRamp2(
        breaks = c(seq(-20,-0.1,length = 19),0,seq(0.1,20,length = 20)),
        colors = c(
            colorRampPalette(c("darkblue","white"))(20)[1:19],
            colorRampPalette(c("white","darkred"))(21)
        )
    )
    
    genecog_ht <- Heatmap(
        name = "Enrichment",
        as.matrix(enrichdf),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col = col_cogs,
        row_names_side = "left",
        column_names_side = "top",
        cell_fun = function(j,i,x,y,width,height,fill){
            if(as.matrix(pvaldf)[i,j] < fisher_pval)
            grid.text(
                "*", x, y, gp = gpar(fontsize = 15)
                ) # add asterisk if p < .05
            }
    )

    if( plotting == TRUE ){
    #    barplot(cogdata_table$freq, col = col_cogs)
    #    pie(cogdata_table$freq, col = col_cogs)
    }

    # Return the results
    res <- list(
        # add here a table with pairmodule-gene-gfam table, long format
        COGperModule = quantdf,
        enrichment = enrichdf,
        pvalue = pvaldf,
        heatmap = genecog_ht
    )

    return(res)
}


