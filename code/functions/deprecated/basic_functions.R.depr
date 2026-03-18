### comparABle
### Basic calculation functions
#' jsd: perform jensen-shannon distance; MUST FIND BETTER WAY
jsd <- function(a,b){
  
  a <- apply(a+1,2,function(x){x/sum(x)})
  b <- apply(b+1,2,function(x){x/sum(x)})
  
  js <- function(p,q){ #jensen-shannon from https://stackoverflow.com/questions/11226627/jensen-shannon-divergence-in-r
    n <- 0.5 * (p + q)
    JS <- sqrt(0.5 * (sum(p * log((p / n)) + sum(q * (log((q / n)))))))
    return(JS)
  }
  
  mat <- matrix(NA,nrow=ncol(a),ncol=ncol(b))
  for (i in 1:ncol(a)){
    for (j in 1:ncol(b)){
      mat[i,j] <- js(a[,i],b[,j])
    }
  }
  
  rownames(mat) <- colnames(a)
  colnames(mat) <- colnames(b)
  
  return(mat)
}

#' topgenes_lmfit
topgenes_lmfit <- function(data, cols_of_interest, top = .1) {
  require(limma)
  require(edgeR)
  
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
  top_genes <- top_genes[top_genes$P.Value < 0.1, ]
  top_genes <- top_genes[order(top_genes$P.Val,top_genes$logFC), ]
  
  # top X genes by FC
  q <- 1 - top
  res <- rownames(top_genes[abs(top_genes$logFC) > quantile(top_genes$logFC,q),])
  
  return(res)
}