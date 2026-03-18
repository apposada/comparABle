#' cooccurrence: implementation of AnaMariaElek/sebepedroslab function, from Levy et al 2023 
#' original function name: `treeFromEnsembleClustering()`
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
cooccurrence <- function(
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
