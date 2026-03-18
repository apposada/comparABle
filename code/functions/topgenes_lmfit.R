#' topgenes_lmfit
topgenes_lmfit <- function(data, cols_of_interest, top = .1, filter_by = "P.Value", full_output = FALSE, only_positive = TRUE) {
  require(limma)
  require(edgeR)
  
  stopifnot(filter_by %in% c("adj.P.Val", "P.Value"))

  # create sample table
  samples = data.frame(
    sample = colnames(data),
    cond = ifelse(colnames(data) %in% cols_of_interest, "X","ctl"),
    row.names = colnames(data)
  )
  
  # create dge object
  dge <- DGEList(
    counts=data,
    samples=samples,
    group=samples$cond
  )
  
  # create design matrix
  design <- model.matrix( ~ samples$cond, dge)
  
  # parse colnames
  colnames(design) <- gsub("samples[$]cond", "", colnames(design))
  colnames(design) <- gsub("samples[$]sample", "", colnames(design))
  names(attr(design,"contrasts")) <- c("condition")
  rownames(design) <- rownames(samples)
  
  # voom fitting
  dge <- calcNormFactors(dge, method="TMM")
  v <- voom(dge, design)
  
  # limma lm fitting
  fitv <- lmFit(v, design)
  fitv <- eBayes(fitv)
  
  # identify significantly upregulated genes
  top_genes <- topTable(fitv, coef="X", n=nrow(data), sort.by="none", adjust.method="BH")
  top_genes_all <- top_genes
  top_genes <- top_genes[top_genes[[filter_by]] < 0.1, ]
  top_genes <- top_genes[order(top_genes[[filter_by]],top_genes$logFC), ]
  top_genes_all

  # top X genes by FC
  q <- 1 - top

  y = rownames(top_genes[abs(top_genes$logFC) > quantile(top_genes$logFC,q),])
  if(only_positive){
    y = rownames(top_genes[top_genes$logFC > 0 & top_genes$logFC > quantile(top_genes$logFC,q),])
  }

  res <-
    list(
      toptable_all = top_genes_all,
      toptable = top_genes,
      topgenes = y
    )

  if(!full_output) res <- y

  return(res)
}
