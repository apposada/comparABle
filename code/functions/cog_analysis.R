### comparABle
### COG analysis

cog_analysis <- function(
    geneset, cogdata,
    col_cogs = c(rainbow(length(LETTERS))[1:25],"black"),
    plotting = TRUE
    ) {
    # Define genes in category; basically subset COG file to retrieve genes present in your dataset of interest
    cogdata_set <- cogdata[cogdata$id %in% geneset, ]

    # Quantify this: make a table()
    cogdata_table <- data.frame(
        table(cogdata_set$cog)
    )
    colnames(cogdata_table) <- c("cog","freq")

    # Add back the COGs that are zero
    cogdata_table <- rbind(
        cogdata_table,
        data.frame(
            cog = LETTERS[which(!(LETTERS %in% cogdata_table$cog))],
            freq = 0
        )
    )
    cogdata_table <- cogdata_table[order(cogdata_table$cog),]

    # Normalise to 1(one)
    cogdata_table_norm <- data.frame(
        cog = cogdata_table$cog,
        freq_norm = cogdata_table$freq / sum(cogdata_table$freq)
    )
    
    if( plotting = TRUE ){
        barplot(cogdata_table$freq, col = col_cogs)
        pie(cogdata_table$freq, col = col_cogs)
    }

    # Return these results
    res <- list(
        genes_cog = cogdata_set,
        cogdata_table = cogdata_table,
        cogdata_table_norm = cogdata_table_norm
    )
    
    return(cogdata_set)

}

#' function to calculate cog enrichment for list  of N group of genes x J cogs.
#' 
cog_enrichment <- function(
    x_modules, x_cog, specific_cogs = FALSE, fisher_pval = 0.05
    ) {

    modules <- sort(unique(x_modules$module))
    if (specific_cogs ! = FALSE) {
        cogs <- sort(unique(specific_cogs)) # must have numerics at the beginning
    } else {
        cogs <- sort(unique(x_cog$cog)) # must have numerics at the beginning
    }

    x_module_cog <- sort(unique(
        merge(
            x_modules,
            x_cog,
            by = 1
        )
    ))

    pvalmat <- matrix(
        rep(1,times = length(modules)*length(cogs)),
        nrow = length(modules),
        ncol = length(cogs)
    )

    enrichmat <- matrix(
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
                    sum(x_module_cog$module == i & x_module_cog$cog ! = j),
                    sum(x_module_cog$module ! = i & x_module_cog$cog == j),
                    sum(x_module_cog$module ! = i & x_module_cog$cog ! = j)
                    ),
                nrow = 2
            )

            enrichmat[h,n] <- 100 * # percent
            (
                contingency[1,1]/(contingency[2,1] + contingency[1,1]) -
                contingency[1,2]/(contingency[2,2] + contingency[1,2])
            )

            pvalmat[h,n] <- fisher.test(contingency)$p.value
        }
    }

    enrichdf <- as.data.frame(enrichmat)
    colnames(enrichdf) <- cogs; rownames(enrichdf) <-  modules

    pvaldf <- as.data.frame(pvalmat)
    colnames(pvaldf) <- cogs ; rownames(pvaldf) <- modules
    pvaldf[pvaldf > fisher_pval ] <- 1 # clipped

    genecog_ht <- Heatmap(
        name = "Gene cog\nEnrichment",
        as.matrix(enrichdf),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col = colorRampPalette(c("#440154","white","#21908C"))(20),
        row_names_side = "left",
        column_names_side = "top",
        cell_fun = function(j,i,x,y,width,height,fill){
            if(as.matrix(pvaldf)[i,j] < fisher_pval)
            grid.text("*", x, y, gp = gpar(fontsize = 15)) # add asterisk if p < .05
            }
    )

    res <- list(
        enrichment = enrichdf,
        pvalue = pvaldf,
        heatmap = genecog_ht
    )

    return(res)
}


