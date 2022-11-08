### comparABle
#' Functions Comparisons across species
#' 
#' Core (essential minimum) input data:
#' - a: expression matrix of species A, quantile-normalized
#' - b: Expression matrix of species B, quantile-normalized
#' - f: Gene Family translation layer between species (GFs)
#' - o: 1:1 orthology translation layer between spp (O)
#' 
#' Basic (extended) input data:
#' - ma: gene modules of species A
#' - mb: gene modules of species B
#' - ga: gene age of spp A
#' - gb: gene age of spp B
#' 
#' Additional parameters:
#' + co: method of correlation (pearson, spearman, JSD)
#' + test: statistical test for module comparison (binom, hypg)
#' + p: subsampling for bootstrapping
#' + n: no. of interations (100)
#' 
#' Basic functions:
#' 
#' To see tidyup and basic calculation functions, check tidyup_functions.R
#' and basic_functions.R
#' 
#' mergedata : merge expression data using 1:1 orthologs
#' mergedata(a, b, o) --> a_o, b_o, ab_o
mergedata <- function(a, b, o) {
  o_a <- o[, c("one2one", "a")]
  o_b <- o[, c("one2one", "b")]
  merge_x_o <- function(x, o) {
    x_o <- merge(x, o, by.x = 0, by.y = 2)
    rownames(x_o) <- x_o$one2one
    x_o <- x_o[, !sapply(x_o, is.character)]
    x_o <- x_o[order(rownames(x_o)), ]
    return(x_o)
  }

  a_o <- merge_x_o(a, o_a)
  b_o <- merge_x_o(b, o_b)
  ab_o <- merge(a_o, b_o, by = 0)
  rownames(ab_o) <- ab_o$Row.names
  ab_o <- ab_o[order(rownames(ab_o)), ]
  ab_o <- ab_o[, !sapply(ab_o, is.character)]

  a_o <- a_o[rownames(a_o) %in% rownames(ab_o), ]
  b_o <- b_o[rownames(b_o) %in% rownames(ab_o), ]

  res <- list(
    a_o = a_o, 
    b_o = b_o, 
    ab_o = ab_o
  )
  return(res)
}
#' 
#' fam_modules: assign gfam to the input modules
#' fam_modules(ma, f) --> df(gene, module, gfam)
#' 
#' 
#' Main functions:
#' 
#' rawcorsp : calculate correlation/distance between spp 
#' using several or user-defined metrics (default: pearson, 
#' spearman, jensen-shannon distance)
#' rawcorsp(a_o, b_o) --> list(pe:matrix, sp:matrix, js:matrix)
rawcorsp <- function(a_o, b_o){ # make one for TFs or genes of interest
  pe <- cor(a_o, b_o, m = "pe")
  sp <- cor(a_o, b_o, m = "sp")
  # je <- jsd(a_o, b_o)
  
  res <- list(
    pe = pe, 
    sp = sp#, 
    # je = jsd
  )
  return(res)
}

#' treecluster: implementation of AnaMariaElek's function, 
#' ALL CREDIT IS HERS. If we ever bundle this as a package, 
#' this should be left out or discussed w her and Arnau...
#' 
#' treecluster(ab_o, ...) --> list(matrix, dendrogram, genes responsible)
#' Build tree by ensemble clustering using downsampling variable genes
#'
#' @param x matrix, genes in rows, cell types (or metacells) in columns
#' @param n numer of iterations of clutering that will be performed 
#'   to calculate coocurence, #commentbyAPP: ideally 1000
#' @param k numeric, number of clusters into which to cut trees 
#'   to get clusters for which the co-occurrence will be calculated
#' @param h numeric, height at which to cut trees to get clusters for
#'   which the co-occurrence will be calculated
#' @param vargenes character, optinal variable genes for sampling
#' @param p numeric, between 0 and 1, fraction of genes for sampling
#' @param clustering_method character, method to use for clustering, 
#'   for available options see method argument in `?hclust()`
#' @param cor_method character, see method argument in `?cor()`
#'

treeFromEnsembleClustering <- function(
  x, n = 1000, k = NULL, h = NULL, vargenes = NULL, p = 0.75,
  bootstrap = FALSE, clustering_algorithm = "hclust",
  clustering_method = "average", cor_method = "pearson"
  ) {
  require(ape)
  ensemble_hc = vector(mode = "list", length = n)
  hs <- rep(h, n)[1:n]
  if (clustering_algorithm == "hclust") {
    for (i in 1:n) { 
      if (bootstrap == TRUE) p = 1
      ids = sample(vargenes, p*length(vargenes), replace = bootstrap)
      hc = hclust(as.dist(1-cor(
        x[ids, ], method = cor_method
      )), clustering_method)
      hc$height <- round(hc$height, 6)
      ensemble_hc[[i]] = cutree(hc, h = hs[i])
    }
  } else if (clustering_algorithm == "nj") {
    for (i in 1:n) {
      if (bootstrap == TRUE) p = 1
      ids = sample(vargenes, p*length(vargenes), replace = bootstrap)
      hc = as.hclust(force.ultrametric(root(nj(as.dist(1-cor(
        x[ids, ], method = cor_method
      ))), outgroup = 1), method = "extend"))
      ensemble_hc[[i]] = cutree(hc, h = hs[i])
    }
  }
  coocmat = matrix(0, nrow = ncol(x), ncol = ncol(x))
  colnames(coocmat) = colnames(x)
  rownames(coocmat) = colnames(x)
  for (i in 1:ncol(x)) {
    for (j in 1:ncol(x)) {
      if (!i<j) {
        a = colnames(x)[i]; b = colnames(x)[j]
        cooc = sum(unlist(lapply(ensemble_hc, function(el) el[a] == el[b])))
        coocmat[i, j] = cooc
        coocmat[j, i] = cooc
      }
    }
  }
  if (clustering_algorithm == "hclust") {
    tree = as.phylo(hclust(dist(coocmat), method = clustering_method))
  } else if (clustering_algorithm == "nj") {
    tree = as.phylo(nj(dist(coocmat)))
  }
  return(list(
    tree = tree, 
    cooccurrence = coocmat
  ))
}

#' commongenes_cor: function to retrieve highly correlating genes
#' across equal-length periods of stages. Far from optimal;
#' workaround until we manage to implement treegapgenes function
commongenes_cor <- function(
  ab_o, across_a = characer(), across_b = character(), 
  min_cor = 0.9) {
  x <- ab_o[ , colnames(ab_o) %in% c(across_a, across_b) ]
  res <- data.frame(
    id = rownames(x), 
    cor = 0
  )
  
  cols_a <- which(colnames(x) %in% across_a)
  cols_b <- which(colnames(x) %in% across_b)
  
  for(i in 1:nrow(x)){
    res$cor[i] <- cor(
      c(as.matrix(x[i, cols_a])),
      c(as.matrix(x[i, cols_b]))
    )
  }
  res <- res[res$cor > min_cor, ]
  #add option to include the original names in spp A and B
  return(res)
}

#' comparemodules: gene family enrichment between modules, 
#' using the gfam object and the modules from each species, 
#' using a binomial statistical test and/or a hypergeometric test.
#' comparemodules(ma, mb, f, ga, gb) --> list(a_f, b_f, ma_f, mb_f, 
#' matrix_hypgeom, matrix_binomial, age_common, age_exclusive)
#' 
#' @param ma association file genes sp a -- gene modules
#' @param mb association file genes sp b -- gene modules
#' @param f association file for gene spp a,b -- gene family
#' 
comparemodules <- function(ma, mb, f){
  fam_of = function(x) {f$gfam[f$id == x]} # lookup function to grab the family of a given gene
  gimme_fams = function(x) { # x is a set of gene names. is this function slow?
    l = c()
    for (i in x){
      l = c(l, fam_of(i))
    }
    l = unique(l)
    return(l)
  }
  gene_POP <- length(unique(f$gfam)) # population size
  
  #' Create a matrix to store stats
  fon = data.frame(
    module_a = "none", 
    module_b = "none", 
    success_in_samples = 0, 
    sample_size = 0, 
    success_in_pop = 0, 
    gene_POP = gene_POP, 
    hypgeom_pval = numeric(1), 
    hypgeom_log = numeric(1), 
    binom_pval =  numeric(1), 
    gfams_common = "none",
    gfams_excl_a = "none",
    gfams_excl_b = "none"
  )
  
  #' Transform into lists for practicality
  ma_list <- split(ma$id, ma$module) 
  mb_list <- split(mb$id, mb$module)
  
  #' Create a matrix to store results (-logpval or similar)
  PVS <- data.frame() # create matrix to store result pval
  PVS_binom <- data.frame()

  
  #' Compute the upper tail of the hypergeometric 
  #' distribution (survival function) for each pair of modules
  #'  as a metric of enrichment #comment from panos @ skarmetalab
  for (i in 1:length(ma_list)){ #need to dramatically optimise speed
    
    a_modulei_name <- names(ma_list[i])
    a_modulei_genes <- ma_list[[i]]
    a_fams <- gimme_fams(a_modulei_genes) #is this slow?
    
    for (j in 1:length(mb_list)){
      
      b_modulej_name <- names(mb_list[j])
      b_modulej_genes <- mb_list[[j]]
      b_fams <- gimme_fams(b_modulej_genes) #is this slow?
      
      sample_size <- length(a_fams)
      success_in_pop <- length(b_fams)
      common_fams <- paste(
        a_fams[which(a_fams %in% b_fams)], collapse = ", "
        )
      exclusive_fams_a <- paste(
        a_fams[which(!(a_fams %in% b_fams))], collapse = ", "
        )
      exclusive_fams_b <- paste(
        b_fams[which(!(b_fams %in% a_fams))], collapse = ", "
        )
      success_in_sample <- length(which(a_fams %in% b_fams))
      
      if (success_in_sample > 0) {
        hypg <- phyper( # HYPGEOMTEST
          # from stackoverflow/questions/8382806/hypergeometric-test-phyper
          q = success_in_sample - 1, # no. of success balls drawn from urn
          m = success_in_pop, # no. of success balls in the urn
          n = gene_POP-success_in_pop, # no. of non-success in the urn
          k = sample_size, # no. of balls drawn from urn
          lower.tail = FALSE
        )
        binom <- binom.test(
          x = success_in_sample, 
          n = sample_size, 
          p = success_in_pop / gene_POP, 
          )$p.value
        
        fon <- rbind(
          fon, 
          c(
            a_modulei_name, b_modulej_name, 
            success_in_sample, sample_size, 
            success_in_pop, gene_POP, as.numeric(hypg), 
            -log(hypg), as.numeric(binom), common_fams, # add here which are the names of the gfams enriched.
            exclusive_fams_a, # exclusive fams in a
            exclusive_fams_b # exclusive fams in b
            )
        )
      } else {
        hypg <- 1
      }
      
      PVS[i, j] <- hypg
      PVS_binom[i, j] <- binom
      
    }
  }
  
  # Tidy up of data 
  rownames(PVS) <- names(ma_list)
  colnames(PVS) <- names(mb_list)
  rownames(PVS_binom) <- names(ma_list)
  colnames(PVS_binom) <- names(mb_list)
  fon <- fon[!(fon$module_a == "none"), ]
  fon$hypgeom_pval <- as.numeric(fon$hypgeom_pval)
  fon$hypgeom_log <- as.numeric(fon$hypgeom_log)
  fon$binom_pval <- as.numeric(fon$binom_pval)
  
  loghypg <- -log(PVS)
  loghypg[loghypg > 30] <- 30
  
  logbinom <- -log(PVS_binom)
  logbinom[logbinom > 30] <- 30
  
  
  res <- list(
    stats = fon, 
    hypgeom = PVS, 
    binom = PVS_binom, 
    loghypg = loghypg, 
    logbinom = logbinom
  )
  
  return(res)
}


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
  gene2go_a, gene2go_b, universe_a, universe_b, universe, 
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
      gene2go_a = gene2go_a, cog_a = cog_a, ma = ma, mb = mb )
  }

  # Analysis of EXCLUSIVE fams
  if (exclusive == TRUE ) {
    keygenes_exclusivefams <- key_genes_in_exclusive_fams(
      stats = x, ga = ga, gb = gb, universe_a = universe_a,
      universe_b = universe_b, gene2go_a = gene2go_a,
      gene2go_b = gene2go_b, cog_a = cog_a, b_cogs = b_cogs,
      ma = ma , mb = mb )
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