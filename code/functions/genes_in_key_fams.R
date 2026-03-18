#' Second part of comparemodules.
#' genes_in_key_fams: Function to retrieve and visualize the enriched
#' gene families by doing gene age enrichment/visualisation;
#' basically it takes the fon table, grab the most significant
#' pairs of comparisons, grab the common gene families and 
#' barplots/pieplots the fraction of genes per gene age, and
#' also grabs the modules from that comparison and performs
#' a gene age enrichment. Perhaps also a gene ontology to
#' see the roles of such families. Or COG.
#' 
#' @param stats data frame output from compare_modules, min format:
#' module_a -- module_b -- pvalue hypgeom -- common gfams -- exclusive gfams
#' @param f association file for gene spp a,b -- gene family
#' @param ma association file genes sp a -- gene modules
#' @param mb association file genes sp b -- gene modules
#' @param top_comparisons number of comparisons to analyse
#' @param common boolean indicating whether to analyse common gfams or not
#' @param exclusive boolean indicating whether to analyse exclusive gfams or not
#' @param same_species boolean indicating if it is intramodular comparisons
#' @param age association file genes sp a -- gene modules
#' @param age_a association file genes sp a -- gene age
#' @param cog_a association file genes sp a -- COG category
#' @param gene2go_a geneID2GO annotation file for GO enrichment (topGO) in sp a
#' @param age_b association file genes sp b -- gene age
#' @param cog_b association file genes sp b -- COG category
#' @param gene2go_b geneID2GO annotation file for GO enrichment (topGO) in sp b
#' @param universe_a gene universe for species a in GO enrichment (topGO)
#' @param universe_b gene universe for species a in GO enrichment (topGO)
#' 
#' 
genes_in_key_fams <- function(
  stats, f, age_a, ma, mb, module_a = FALSE, module_b = FALSE, # add which metric: hypgeom or binomial
  top_comparisons = 10, common = TRUE, exclusive = FALSE , 
  age_b, cog_a, cog_b, same_species = FALSE, 
  gene2go_a, gene2go_b, universe_a, universe_b, universe, sep,
  ...
){
  #' The proper way to do this is, taking the exact genes
  #' from that module based on what gene families they
  #' are from.
  #' 
  #' If gene123 from the gfam ABCD is found significantly
  #' present in the two species, you do not want to take all
  #' the genes from spA that belong to the gfam ABCD. You
  #' want to take gene123 AND any other gene from ABCD also
  #' expressed in the same module as gene123. For that, 
  #' you need to integrate information from module
  #' membership, gene age, and gfam.
  #' 
  #' You basically want to revert to the geneid_moduleid
  #' type of information but JUST with the genes from
  #' families found important.
  #' 
  #' How to do this?
  #' Take common gfam, take all genes from those fams, 
  #' subset for genes found in module of interest module_a, 
  #' --> thus obtaining the list of genes from module_a and 
  #' common gfams. THEN you can do GO analysis, COG, and
  #' gene age...
  
  stats$hypgeom_log <- as.numeric(stats$hypgeom_log)
  
  if( (module_a !=  F) & (module_b !=  F)){
    x <- stats[
      stats$module_a == module_a & stats$module_b == module_b, 
    ]
  } else{
    x <- stats[
      rev(order(stats$hypgeom_log)), 
    ]
    if( same_species != FALSE ) {
      x <- x[x$module_a != x$module_b,] # this is for plei, same species, same modules... this needs to be implemented somehow
    }
    x <- x[1:top_comparisons, ]
  }
  
  # Analysis of COMMON fams
  if (common == TRUE ) {
    keygenes_commonfams <- key_genes_in_common_fams(
      stats = x, ga = ga, universe = universe,
      gene2go_a = gene2go_a, cog_a = cog_a, ma = ma, mb = mb,
      sep = sep)
  }
  
  # Analysis of EXCLUSIVE fams
  if (exclusive == TRUE ) {
    keygenes_exclusivefams <- key_genes_in_exclusive_fams(
      stats = x, ga = ga, gb = gb, universe_a = universe_a,
      universe_b = universe_b, gene2go_a = gene2go_a,
      gene2go_b = gene2go_b, cog_a = cog_a, b_cogs = b_cogs,
      ma = ma , mb = mb, sep = sep)
  }
  
  # Organise results
  
  if( common == TRUE ){
    
    if (exclusive == TRUE) {
      
      res <- list(
        commonfams = keygenes_commonfams,
        exclusivefams = keygenes_exclusivefams
      )
      
    }
    
    res <- list(
      commonfams = keygenes_commonfams
    )
    
  }
  
  # return outputs
  return(res)
}

#' js_with_subsampling: applies bootstrap using the j-s function
#' 
jsd_with_subsampling <- function(
  a_o,b_o, n = 1000, p = 0.75, bootstrap = FALSE
) {
  
  ensemble_js = vector(mode = "list", length = n)
  
  genes = rownames(a_o)[rownames(a_o) %in% rownames(b_o)]
  
  for (i in 1:n) { 
    if (bootstrap == TRUE) p = 1
    
    ids = sample(genes, p*length(genes), replace = bootstrap)
    
    a_ <- a_o[ids,]
    b_ <- b_o[ids,]
    
    js <- jsd(a_, b_)
    
    ensemble_js[[i]] = js
  }
  
  # calculate the mean matrix
  mean_js <- Reduce("+", ensemble_js) / n
  
  # sd_js <- Reduce("+", lapply(ensemble_js, function(x) (x - apply(ensemble_js, 2, mean))^2)) / n-1
  # sd_js <- sqrt(sd_js)
  
  rownames(mean_js) = colnames(a_o)
  colnames(mean_js) = colnames(b_o)
  
  # rownames(sd_js) = colnames(a_o)
  # colnames(sd_js) = colnames(b_o)
  
  res <- list(
    mean = mean_js#,
    # sd = sd_js
  )
  
  return(res)
  
}
