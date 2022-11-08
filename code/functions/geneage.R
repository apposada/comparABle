### gene age enrichments
#' function to calculate gene age enrichment for list  of N group of genes x J 
#' gene ages or phylostrata. From Andrej Ondracka and A.Perez-Posada ca. 2018
#'

require(ComplexHeatmap)
require(circlize)
require(ggplot2)
require(dplyr)
#' 
#' @param x_modules association file between genes -- modules/categories/clusters...
#' @param x_cog association file between genes -- COG categories (any category)
#' @param specific_cogs either FALSE or vector of specific COGs to be tested
#' @param fisher_pval minimum pvalue to consider significant in the heatmap visualisation
#' @param potting whether to plot the heatmap or not
#' @param functional_categories either FALSE or a table with alternative names for COG categories
#' 

gene_age_enrichment <- function(
    x_modules, x_age, phylostrata = FALSE, fisher_pval = 0.05
    ) {

    modules <- sort(unique(x_modules$module))
    if (phylostrata != FALSE) {
        ages <- sort(unique(phylostrata)) # must have numerics at the beginning
    } else {
        ages <- sort(unique(x_age$age)) # must have numerics at the beginning
    }

    x_module_age <- unique(
        merge(
            x_modules,
            x_age,
            by = 1
        )
    )

    quantmat <- matrix(
        rep(0,times = length(modules)*length(ages)),
        nrow = length(modules),
        ncol = length(ages)
    )

    enrichmat <- matrix(
        rep(1,times = length(modules)*length(ages)),
        nrow = length(modules),
        ncol = length(ages)
    )

    pvalmat <- matrix(
        rep(1,times = length(modules)*length(ages)),
        nrow = length(modules),
        ncol = length(ages)
    )

    for (i in modules) {
        h <- which(modules == i)
        for (j in ages) {
            n <- which(ages == j)

            contingency <- matrix(
                c(
                    sum(x_module_age$module == i & x_module_age$age == j),
                    sum(x_module_age$module == i & x_module_age$age != j),
                    sum(x_module_age$module != i & x_module_age$age == j),
                    sum(x_module_age$module != i & x_module_age$age != j)
                    ),
                nrow = 2
            )

            quantmat[h,n] <- contingency[1,1] # no. of COGj-genes in module i

            enrichmat[h,n] <- 100 * # percent
            (
                contingency[1,1]/(contingency[2,1] + contingency[1,1]) -
                contingency[1,2]/(contingency[2,2] + contingency[1,2])
            )

            pvalmat[h,n] <- fisher.test(contingency)$p.value
        }
    }

    quantdf <- as.data.frame(quantmat)
    colnames(quantdf) <- ages; rownames(quantdf) <- modules

    enrichdf <- as.data.frame(enrichmat)
    colnames(enrichdf) <- ages; rownames(enrichdf) <- modules

    pvaldf <- as.data.frame(pvalmat)
    colnames(pvaldf) <- ages ; rownames(pvaldf) <- modules
    pvaldf[pvaldf > fisher_pval ] <- 1 # clipped

    col_age <- colorRamp2(
        breaks = c(seq(-20,-0.1,length = 19),0,seq(0.1,20,length = 20)),
        colors = c(
            colorRampPalette(c("#440154","white"))(20)[1:19],
            colorRampPalette(c("white","#21908C"))(21)
            )
        )
    geneage_ht <- Heatmap(
        name = "Gene age\nEnrichment",
        as.matrix(enrichdf),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col = col_age,
        row_names_side = "left",
        column_names_side = "top",
        cell_fun = function(j,i,x,y,width,height,fill){
            if(as.matrix(pvaldf)[i,j] < fisher_pval)
            grid.text(
                "*", x, y, gp = gpar(fontsize = 15)
                ) # add asterisk if p < .05
            }
    )

    res <- list(
        # add here a pairmodule-gene-gfam table, longformat
        AgeperModule = quantdf, # this can be normalised and it makes for 
        enrichment = enrichdf,
        pvalue = pvaldf,
        heatmap = geneage_ht
    )

    return(res)
}


#' 
#' It would be useful to have a phylostrata dataset not only
#' for species-specific genes but also for gene families;
#' eg:
#' gfam_001234  Mzoa
#' gfam_001235  Opisthokonta
#' gfam_001236  Bilateria
#' ...
#' 
#' see if this can be retrieved from OF/OMA/whatever?
#' This would make everything SO SO MUCH EASIER. Could be
#' easily interchangeable with genes/gfam levels
#' 
#' 
#' Import equivalence of OMA_ and integer fams
#' 
#' 
#' Import an orthogroups file to calculate gene age using dollo parsimony
#' 
#' 
#' Import a dollo parsimony table, assign the gene ages
#' 
#' 
#' Assign gene ages to all the genes in the datasets of spp a and spp b
#' 
#' 
#' 
# c("#f1d1b5", "#f8b195", "#f67280", "#c06c84", "#6c5b7b", "#568ea6", "#305f72", "#355c7d", "#F9A26C", "#f67280", "#801336", "#6c5b7b", "#2F9599", "#355C7D")