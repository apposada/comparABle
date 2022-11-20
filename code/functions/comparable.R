comparABle <- function(
  a_name,
  b_name,
  a,
  b,
  o,
  f,
  ma,
  mb,
  ga,
  gb,
  cog_a,
  cog_b,
  across_a,
  across_b,
  cooc_n = 1000,
  cooc_h = c(0.75,0.95),
  cooc_clustering_algorithm = "hclust",
  cooc_clustering_method = "average",
  cooc_cor_method = "pearson",
  cooc_p = 0.1,
  cooc_vargenes = rownames(merge_ab$ab_o),
  highlyvariable = TRUE,
  ...
  ) {
  
  # TidyUp
  print("Tidy up data")
  samples_a = levels(condition_x)
  if (any(apply(a,2,is.character)) == TRUE) {
    a = a[,!(sapply(a, is.character))]
  }
  a = qnorm(a)
  a = tidyup(a, highlyvariable = highlyvariable) # remove genes with 0 tpms
  a = rep2means(samples_a,a)
  
  if (any(apply(b,2,is.character)) == TRUE) {
    b = b[,!(sapply(b, is.character))]
  }
  samples_b = unique(sub("_.$", "", colnames(b)))
  b = qnorm(b)
  b = tidyup(b, highlyvariable = highlyvariable)
  b = rep2means(samples_b,b) # remove genes with 0 tpms
  
  o = pair2id(o)
  
  colnames(ma) <- c("id","module")
  colnames(mb) <- c("id","module")
  
  # MERGE
  print ("Merge data")
  merge_ab <- mergedata(a,b,o)
  
  # CORRELATIONS
  print("Correlations")
  cors <- rawcorsp(merge_ab$a_o,merge_ab$b_o) # FIX JSD
  
  print("PCA")
  pi <- prcomp(t(merge_ab$ab_o))
  
  # CO-OCCURENCE MATRIX
  print("Co-Occurrence")
  set.seed(4343)
  cooc <- treeFromEnsembleClustering(
    x=merge_ab$ab_o, p=cooc_p, h=cooc_h,  n = cooc_n, vargenes = rownames(merge_ab$ab_o), bootstrap=FALSE,
    clustering_algorithm=cooc_clustering_algorithm, clustering_method=cooc_clustering_method, 
    cor_method=cor_method
  )
  
  # Co-occurrence heatmap
  cooc_hm <- Heatmap(
    cooc$cooccurrence,
    col = colorRamp2(
      c(seq(min(cooc$cooccurrence),
            max(cooc$cooccurrence),
            length=9
      )
      ),
      colors=c(
        c("white",'#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#990000')
      )
    ), 
    name="co-occurence"
  )
  
  # COMMON GENES IN CORRELATIONS
  print("Common genes in Correlations")
  # This should be determined autommatically ... or with an "if (stages_hicor are known) {} else {} ..."
  ab_common_genes_cor <- commongenes_cor(
    merge_ab$ab_o,
    across_a = across_a, 
    across_b = across_b,
    min_cor = 0.8
    )
  
  ab_common_genes_cor <- merge(
    ab_common_genes_cor,
    o,
    by = 1
  )
  
  print("Common genes in Correlations (GO)")
  ab_common_genes_cor_GOs <- getGOs(
    genelist = list(a = ab_common_genes_cor$a),
    gene_universe = pfla_all_gene_names,
    gene2GO = pfla_geneID2GO
    )
  
  print("Common genes in Correlations (age)")
  ab_common_genes_cor_age <- gene_age_enrichment(
    x_modules = data.frame(
      id = pfla_all_gene_names,
      module = ifelse(
        pfla_all_gene_names %in% ab_common_genes_cor$a,
        "common",
        "non_common"
      )
    ),
    x_age = ga[!(ga$age %in% c("Ambu","Hemich","Pfla_specif")),] # custom vector should be a variable of common evol.nodes
  )
  
  # COMPARE MODULES
  print("Pairwise Orthology Overlap Strategy across modules -- hypergeometric and binonmial tests")
  modulecomp_ab <- comparemodules(ma,mb,f)
  
  # Genes in key fams across modules
  print("Getting info on genes from shared families across modules")
  ab_common_genes_details <- genes_in_key_fams(
    stats = modulecomp_ab$stats,
    top_comparisons = 20,
    f = f,
    ma = ma,
    mb = mb,
    age_a = ga[!(ga$age %in% c("Ambu","Hemich","Pfla_specif")),], # custom vector should be a variable of common evol.nodes
    cog_a = cog_a,
    gene2go_a = pfla_geneID2GO,
    universe_a = pfla_all_gene_names,
    universe = pfla_all_gene_names,
    sep = ",",
    module_a = FALSE,
    module_b = FALSE,
    common = TRUE,
    exclusive = FALSE , 
    same_species = FALSE
    )
  
  
  # PLOTS
  print("Plots")
  # STORE ALL PLOTS IN FUNCTIONS AS DONE BY @cartwheel ON https://stackoverflow.com/questions/29583849/save-a-plot-in-an-object
  # PCA
  ab_pca <- function(){plot_pca_ab(pca = pi,ab_o = merge_ab)}
  
  # Correlation, Coocurrence
  ab_spearman <- Heatmap(
    column_title = paste0("correlation ",a_name," vs ",b_name),
    cors$sp,
    name = "spearman",
    cluster_rows = F,
    cluster_columns = F,
    col = rev(sequential_hcl(10,"YlOrRd"))
  )
  
  ab_pearson <- Heatmap(
    cors$pe,
    name = "correlation",
    cluster_rows = F,
    cluster_columns = F,
    col = rev(sequential_hcl(10,"YlOrRd"))
  )
  
  ab_jsd <- Heatmap(
    cors$js,
    name = "correlation",
    cluster_rows = F,
    cluster_columns = F,
    col = rev(sequential_hcl(10,"YlOrRd"))
  )
  
  ab_cooc_hm <- function() {
    heatmap(
      cooc$cooccurrence,
      symm = T,
      Rowv = cooc$tree$edge
    )
  }
  
  # Common genes, highly correlated
  high_correlation_genes <- function(){
    par(mfrow = c(1,2))
    boxplot(
      main = a_name,
      ab_common_genes_cor[,3:6],
      col = viridis::inferno(7)[3:7],
      xlab = "stage",
      ylab = "scaled expression"
    )
    boxplot(
      main = a_name,
      ab_common_genes_cor[,7:10],
      col = viridis::magma(7)[3:7],
      xlab = "stage",
      ylab = "scaled expression"
    )
    par(mfrow = c(1,1))
  }
  
  # GO terms
  #' turn this into a GGplot, same for the whole plethora of
  #' enriched GO terms
  # 
  # plot(
  #   x = 
  #     ab_common_genes_GOs$GOtable$a$Significant /
  #     ab_common_genes_GOs$GOtable$a$Expected,
  #   y = 
  #     -log(
  #       as.numeric(
  #         ab_common_genes_GOs$GOtable$a$classicFisher
  #       )
  #     ),
  #   xlab = "FC no. Obs/Exp genes",
  #   ylab = "-logpvalue"
  # )
  # 
  # Age of highly cor genes
  high_correlation_genes_age <- function(){
    par(mfrow = c(1,2))
    barplot(
      t(t(ab_common_genes_cor_age$enrichment[1,c(5,3,4,1,2)]))
    )
    barplot(t(t(ab_common_genes_cor_age$AgeperModule[1,c(5,3,4,1,2)])))
    par(mfrow = c(1,1))
  }
  
  
  # Comparison of modules, binomial test
  # turn this into complexheatmap
  modulecomp_ab_logbinom_hm <- Heatmap(
    modulecomp_ab$logbinom,
    col = rev(sequential_hcl(10,"Blues 3")),
    cluster_columns = F,
    cluster_rows = F
  )
  
  # Comparison of modules, hypergeometric test
  modulecomp_ab_loghypg_hm <- Heatmap(
    modulecomp_ab$loghypg,
    col = rev(sequential_hcl(10,"YlOrBr")),
    cluster_cols = F,
    cluster_rows = F
  )
  
  print("Generating results")
  
  res <- list(
    input = list(
      a,
      b,
      o,
      f,
      ma,
      mb,
      ga,
      cog_a
    ),
    merged_data = merge_ab,
    pairwise_correlations = cors,
    pca_analysis = pi,
    coocurrence_analysis = cooc,
    high_corr_genes = list(
      table = ab_common_genes_cor,
      GOs = ab_common_genes_cor_GOs,
      age = ab_common_genes_cor_age
    ),
    orthology_overlap_modules = list(
      pairwise_module_comparison = modulecomp_ab,
      genes_in_common_fams = ab_common_genes_details
    ),
    plots = list( # these functions will NOT work outside the function call I think...
      pca = ab_pc(),
      spearman_cor = ab_spearman,
      pearson_cor = ab_pearson,
      jensen_shannon = ab_jsd,
      coocurrence_hm = ab_cooc_hm(),
      coocurrence_Heatmap = cooc_hm,
      high_correlated_genes_boxplot = high_correlation_genes(),
      high_correlated_genes_boxplot_age = high_correlation_genes_age(),
      orthology_overlap_binomial_hm = modulecomp_ab_logbinom_hm,
      orthology_overlap_hypgeom_hm = modulecomp_ab_loghypg_hm
    )
  )
  
  return(res)
}