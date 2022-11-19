key_genes_in_common_fams <- function(...){ #IMPORTANT CAVEAT SEE BELOW
    x_fams_a <- list()
    x_comparison_modules <- data.frame( # rethink the name of this
    id = "none",
    module = "none"
    )
    x_fams_b <- list()

    for (i in 1:nrow(x)){ # for each instance in the stats tab; ie for each comparison
    mod_a = x$module_a[i] # module in spp a of that comparison
    fams_i_common <- unlist(strsplit(x$gfams_excl_a[i], ", ")) # common fams of that comparison
    genes_module_a <- ma$id[ma$module == mod_a] # genes of spp a in module a of that comparison
    genes_module_a_commonfams <- genes_module_a[ #those genes
        genes_module_a %in% # but only those in
        f[ # the table of gene fams
            f$gfam %in% fams_i_common,  # whose associated fam is among the common ones
            2 # column 2 for the genes
        ]
    ]
    x_fams_a[[i]] <- genes_module_a_commonfams
    x_comparison_modules <- rbind(
        x_comparison_modules,
        data.frame(
        id = genes_module_a_commonfams , 
        module = paste(
            x$module_a[i],
            x$module_b[i],
            sep = "__"
        )
        )
    )

    mod_b = x$module_b[i]
    fams_i_common <- unlist(strsplit(x$gfams_excl_b[i], ", "))
    genes_module_b <- mb$id[mb$module == mod_b]
    genes_module_b_commonfams <- genes_module_b[
        genes_module_b %in% 
        f[
            f$gfam %in% fams_i_common, 
            2 # column containing the genes
        ]
    ]
    x_fams_b[[i]] <- genes_module_b_commonfams
    }

    names(x_fams_a) <- paste(
    x$module_a,
    x$module_b,
    sep = "__" # this will allow easier parsing
    )

    names(x_fams_b) <- paste(
    x$module_a,
    x$module_b,
    sep = "__" # this will allow easier parsing
    )

    #' ... Alternatively, one might say that this is abt gene 
    #' families now and not the info on a particular species.
    #' Meaning it is necessary a prior analysis of COG and
    #' stuff for EVERY gene family and then subset that
    #' analysis to speak about the particular gene families?

    #gene age bar/pieplot
    #barplot; check which one looks better . maybe a grid of pie charts?
    #maybe a grid of piecharts of ALL the comparisons and only in color/highlighted those that are significant?? does Heatmap() allow this?
    barplot( 
    table(
        x_comparison_modules$module
    ),
    las = 2
    )

    #gene age enrichment (barplot of FC up--down )
    # IMPORTANT: this MUST be reimplemented because at the moment it is performing comparisons of enrichment BETWEEN the genes common across pair of modules, not of every set of common genes against the whole gene set of organism a or against the set of all genes in module x of species a. Perhaps the latter is more informative.
    commonfams_age <- gene_age_enrichment(
        # perhaps this should be replaced for another function that retrieves the gene age enrichment between element i in the list of genes/modules, against a certain set of genes (user defined, or both: all genes of a, and on the other side all genes of a belonging to module of interest)
    x_modules = x_comparison_modules,
    x_age = ga,
    phylostrata = phylostrata
    )

    #heatmap of %orthogroups per gene age of relevance
    #' table of #gfams per module per age and normalise. Make a heatmap. 
    #' Will implement in the future

    #cog_enrichment
    # IMPORTANT: this MUST be reimplemented because at the moment it is performing comparisons of enrichment BETWEEN the genes common across pair of modules, not of every set of common genes against the whole gene set of organism a or against the set of all genes in module x of species a. Perhaps the latter is more informative.
    x_fams_a_COGs <- cog_enrichment_analysis(
    # same as above
    x_modules = x_comparison_modules,
    x_cog = a_cogs,
    specific_cogs = specific_cos
    )

    #topgo
    x_fams_a_GOs <- list()
    for (i in 1:nrow(x_fams_a)){
    x_fams_a_GOs[[i]] <- getGOs(x_fams_a[[i]],gene_universe = universe)
    }

    #Data return

    res_common = list(
        table_a_common = x_comparison_modules,
        list_a_common = x_fams_a,
        age_a_common = commonfams_age,
        cog_a_comon = x_fams_COGs,
        go_a_common = x_fams_GOs
    )

    return(res_common)
}



### Genes in exclusive gene families
key_genes_in_exclusive_fams <- function(...){
    x_fams_a_exclusive <- list() # list format
    x_comparison_modules_a_exclusive <- data.frame( # table format # rethink the name of this
    id = "none",
    module = "none"
    )
    x_fams_b_exclusive <- list()
    x_comparison_modules_b_exclusive <- data.frame( # table format # rethink the name of this
    id = "none",
    module = "none"
    )


    # This down below will probably turn into a basic function
    for (i in 1:nrow(x)){ # for each instance in the stats tab; ie for each comparison
        mod_a = x$module_a[i] # module in spp a of that comparison
        fams_i_a_excl <- unlist(strsplit(x$gfams_common[i], ", ")) # common fams of that comparison
        genes_module_a <- ma$id[ma$module == mod_a] # genes of spp a in module a of that comparison
        genes_module_a_commonfams <- genes_module_a[ #those genes
            genes_module_a %in% # but only those in
            f[ # the table of gene fams
                f$gfam %in% fams_i_a_excl,  # whose associated fam is among the common ones
                2 # column 2 for the genes
            ]
        ]
        x_fams_a_exclusive[[i]] <- genes_module_a_commonfams # list format
        x_comparison_modules_a_exclusive <- rbind( # table format
            x_comparison_modules_a_exclusive,
            data.frame(
            id = genes_module_a_commonfams , 
            module = paste(
                x$module_a[i],
                x$module_b[i],
                sep = "__"
            )
            )

        mod_b = x$module_b[i] # module in spp a of that comparison
        fams_i_b_excl <- unlist(strsplit(x$gfams_common[i], ", ")) # common fams of that comparison
        genes_module_b <- mb$id[mb$module == mod_b] # genes of spp a in module a of that comparison
        genes_module_b_commonfams <- genes_module_b[ #those genes
            genes_module_b %in% # but only those in
            f[ # the table of gene fams
                f$gfam %in% fams_i_b_excl,  # whose associated fam is among the common ones
                2 # column 2 for the genes
            ]
        ]
        x_fams_b_exclusive[[i]] <- genes_module_b_commonfams # list format
        x_comparison_modules_b_exclusive <- rbind( # table format
            x_comparison_modules_b_exclusive,
            data.frame(
            id = genes_module_b_commonfams , 
            module = paste(
                x$module_a[i],
                x$module_b[i],
                sep = "__"
            )
            )

        )
    }

    names(x_fams_a_exclusive) <- paste(
        x$module_a,
        x$module_b,
        sep = "__" # this will allow easier parsing
    )

    names(x_fams_b_exclusive) <- paste(
        x$module_a,
        x$module_b,
        sep = "__" # this will allow easier parsing
    )

    barplot( 
        table(
            x_comparison_modules_a_exclusive$module
        ),
        las = 2
    )

    barplot( 
        table(
            x_comparison_modules_b_exclusive$module
        ),
        las = 2
    )

    #gene age enrichment (barplot of FC up--down )
    commonfams_age_a_exclusive <- gene_age_enrichment(
    x_modules = x_comparison_modules_a_exclusive,
    x_age = ga
    )

    commonfams_age_b_exclusive <- gene_age_enrichment(
    x_modules = x_comparison_modules_b_exclusive,
    x_age = gb
    )


    #heatmap of %orthogroups per gene age of relevance
    #' table of #gfams per module per age and normalise. Make a heatmap. 
    #' Will implement in the future

    #cog_enrichment
    x_fams_a_exclusive_COGs <- cog_enrichment_analysis(
    x_modules = x_comparison_modules_a_exclusive,
    x_cog = a_cogs,
    specific_cogs = specific_cos
    )

    x_fams_b_exclusive_COGs <- cog_enrichment_analysis(
    x_modules = x_comparison_modules_b_exclusive,
    x_cog = b_cogs,
    specific_cogs = specific_cos
    )

    #topgo
    x_fams_a_exclusive_GOs <- list()
    for (i in 1:nrow(x_fams_a_exclusive)){
    x_fams_a_exclusive_GOs[[i]] <- getGOs(x_fams_a_exclusive[[i]],gene_universe = universe_a)
    x_fams_b_exclusive_GOs <- list()
    for (i in 1:nrow(x_fams_b_exclusive)){
    x_fams_b_exclusive_GOs[[i]] <- getGOs(x_fams_b_exclusive[[i]],gene_universe = universe_b)
    }

    res_excl_a = list(
        table_a_exclusive = x_comparison_modules_a_exclusive,
        list_a_exclusive = x_fams_a_exclusive,
        age_a_exclusive = commonfams_age_a_exclusive,
        cog_a_exclusive = x_fams_a_exclusive_COGs,
        go_a_exclusive = x_fams_a_exclusive_GOs
    )

    res_excl_b = list(
        table_b_exclusive = x_comparison_modules_b_exclusive,
        list_b_exclusive = x_fams_b_exclusive,
        age_b_exclusive = commonfams_age_b_exclusive,
        cog_b_exclusive = x_fams_b_exclusive_COGs,
        go_b_exclusive = x_fams_b_exclusive_GOs
    )

    return(
        res_excl_a,
        res_excl_b
    )
}