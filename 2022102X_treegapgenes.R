#' Get gap genes for internal nodes in a tree
#' 
#' @param tree
#' @param feature_matrix, rows should be features (genes or motifs) and columns should be 
#'   tips of the tree; data should be log2 scaled
#' @param feature_inclusion_ths numeric, threshold to consider a feature (gene, motif) (default: 1)
#' @param branch_length_ths numeric, threshold for selecting long branches (default: 0.1)
#' @param feature_in_thrs min threshold median value for feature in 
#'   the columns of feature_matrix in the selected node split
#' @param feature_out_thrs max threshold median value for features in 
#'   the columns of feature_matrix outside of selected node split
#' @param method,methodbg character, "absolute" or "median"
#' @param ncores integer (default: `detectCores()-1`)
#' @param verbose logical (default: FALSE)
#' 

tree_gap_genes <- function(
  tree, feature_matrix, branch_length_thrs=0.1,
  feature_in_thrs = 1.5, feature_out_thrs = 1,
  method="absolute", methodbg="absolute", abs_leakyness=0, abs_leakynessbg=0.05,
  ncores=detectCores()-1, verbose=FALSE
) {
  
  # get tree data
  tree_tb=as_tibble(tree)
  treet=tidytree::as.treedata(tree)
  
  # get the order of the tips
  tip_labels <- treet@phylo$tip.label
  is_tip <- treet@phylo$edge[,2] <= length(tip_labels)
  ordered_tips <- treet@phylo$edge[is_tip, 2]
  
  # get the ordered tip labels
  tip_labels <- treet@phylo$tip.label[ordered_tips]
  
  # select iternal nodes with large branch lengths. CRITICAL!!
  # (start at n_tips + 2 (because the first internal node is the whole tree!)
  nodes <- 1:nrow(tree_tb)
  nodes_above_ths <- tree_tb$branch.length[(length(tip_labels)+2):nrow(tree_tb)] > branch_length_thrs
  long_nodes <- which(nodes_above_ths)+length(tip_labels)+1
  
  calc_nodes=c(ordered_tips,long_nodes)
  
  # matrix to store the binary split info  1/0
  m_splits=matrix(0,ncol=length(tip_labels),nrow=length(calc_nodes))
  colnames(m_splits)=tip_labels
  rownames(m_splits)=calc_nodes #1:nrow(m_splits)
  # list to store the enriched feature in the tips under the node vs all other tips
  top_features=vector("list",length=length(nodes))
  # list to store the enriched features in the tips under the node vs sister clade (branching from mrca)
  top_feature_sister=vector("list",length=length(nodes))
  # list to store the anti-enriched features (or enriched oustide) the node
  top_feature_anti=vector("list",length=length(nodes))
  
  # require(foreach)
  # require(doParallel)
  # maxcores <- detectCores()
  # if (ncores > maxcores) {
  #   warning(sprintf(
  #     "Specified numer of cores (%s) exceeds the number of available cores (%s)!\nContinuing using %s cores.",
  #     ncores, maxcores, maxcores
  #   ))
  #   ncores=maxcores
  # }
  #registerDoParallel(ncores)
  #i=0
  #foreach (node=nodes) %dopar% {
  for(node in nodes){
    if (node %in% calc_nodes) {
      if (verbose==TRUE) print(sprintf("Starting analysis for node %s", node))
      #i=i+1
      i=match(node,calc_nodes)
      node_parent <- treeio::parent(tree,node)
      node_children <- treeio::child(tree,node_parent)
      node_sibling <- setdiff(node_children,node)
      if (node %in% ordered_tips) {
        tips_in <- tip_labels[node]
      } else {
        tips_in <- treeio::tree_subset(treet,node,levels_back=0)@phylo$tip.label
      }
      if (verbose==TRUE) print(sprintf("Tips in: %s", paste(tips_in,collapse=", ")))
      sister_label <- tree_tb[tree_tb$node%in%node_sibling,]$label
      tips_sister <- tryCatch({
        if (any(!is.na(sister_label))) {
          sister_label[!is.na(sister_label)]
        } else {
          unlist(lapply(node_sibling, function(node_sibling) 
            treeio::tree_subset(treet,node_sibling,levels_back=0)@phylo$tip.label
          ))
        }
      })#, error=function(e) NULL)
      if (verbose==TRUE) print(sprintf("Sister tips: %s", paste(tips_sister,collapse=", ")))
      tips_out <- setdiff(treet@phylo$tip.label,tips_in)
      if (verbose==TRUE) print(sprintf("Out tips: %s", paste(tips_out,collapse=", ")))
      
      # fill the reference split matrix
      m_splits[i,tips_in] <- 2
      m_splits[i,tips_out] <- 1
      .calc_leakyness <- function(abs_leakyness,n) {
        if (abs_leakyness<1) {
          pmin(round(c(1-abs_leakyness)*n), n)
        } else {
          pmin(round(n-abs_leakyness), n)
        }
      }
      if (method=="median") {
        fin=apply(feature_matrix[,tips_in,drop=F],1,median) > feature_in_thrs
        fin_inv=apply(feature_matrix[,tips_out,drop=F],1,median) > feature_in_thrs
      } else if (method=="absolute") {
        fin=apply(feature_matrix[,tips_in,drop=F], 1, function(x) 
          #all(x > feature_in_thrs)
          !(sum(x > feature_in_thrs) < .calc_leakyness(abs_leakyness,n=length(tips_in)))
        )
        fin_inv=apply(feature_matrix[,tips_out,drop=F], 1, function(x) 
          #all(x > feature_in_thrs)
          !(sum(x > feature_in_thrs) <  .calc_leakyness(abs_leakyness,n=length(tips_out)))
        )
      }
      if (methodbg=="median") {
        fout=apply(feature_matrix[,tips_out,drop=F],1,median) < feature_out_thrs
        fout_inv=apply(feature_matrix[,tips_in,drop=F],1,median) < feature_out_thrs
        fout_sister=apply(feature_matrix[,tips_sister,drop=F],1,median) < feature_out_thrs
      } else if (methodbg=="absolute") {
        fout=apply(feature_matrix[,tips_out,drop=F], 1, function(x) 
          #all(x < feature_out_thrs)
          !(sum(x < feature_out_thrs) < .calc_leakyness(abs_leakynessbg,n=length(tips_out)))
        )
        fout_inv=apply(feature_matrix[,tips_in,drop=F], 1, function(x) 
          #all(x < feature_out_thrs)
          !(sum(x < feature_out_thrs) < .calc_leakyness(abs_leakynessbg,n=length(tips_in)))
        )
        fout_sister=apply(feature_matrix[,tips_sister,drop=F], 1, function(x) 
          #all(x < feature_out_thrs)
          !(sum(x < feature_out_thrs) < .calc_leakyness(abs_leakynessbg,n=length(tips_sister)))
        )
      }
      
      f_in=which(fin & fout)
      f_sister=which(fin & fout_sister)
      f_out=which(fin_inv & fout_inv)
      
      # now add each feature to the corresponding list (or "none")
      if(length(f_in)>0){
        top_feat=rownames(feature_matrix)[f_in]
        if(length(f_in)>1)
          top_feat=top_feat[order(apply(feature_matrix[top_feat,tips_in,drop=F],1,median),decreasing=T)] 
        top_features[[node]]=top_feat
      } else { top_features[[node]]="NONE"}	
      
      if(length(f_sister)>0){
        top_feat=rownames(feature_matrix)[f_sister]
        if(length(f_sister)>1)
          top_feat=top_feat[order(apply(feature_matrix[top_feat,tips_in,drop=F],1,median),decreasing=T)] 
        top_feature_sister[[node]]=top_feat
      } else { top_feature_sister[[node]]="NONE"}	
      
      if(length(f_out)>0){
        top_anti=rownames(feature_matrix)[f_out]
        if(length(f_out)>1)
          top_anti=top_anti[order(apply(feature_matrix[top_anti,tips_out,drop=F],1,median),decreasing=T)]
        top_feature_anti[[node]]=top_anti
      } else { top_feature_anti[[node]]="NONE"}
      
      if (verbose==TRUE) print(sprintf("Finished analysis for node %s", node))
      
    } else {
      if (verbose==TRUE) 
        print(sprintf(
          "Node %s is shorther than threshold branch length (%s); skipping.", 
          node, branch_length_thrs
        ))
      top_features[[node]] <- "NOT CALCULATED"
      top_feature_sister[[node]] <- "NOT CALCULATED"
      top_feature_anti[[node]] <- "NOT CALCULATED"
    }
  }
  #stopImplicitCluster()
  hm <- Heatmap(
    m_splits, col=c("0"="gray88","1"="lightgrey","2"="black"), show_heatmap_legend=FALSE,
    border = TRUE, rect_gp = gpar(col="gray80"),
    cluster_columns = FALSE, cluster_rows = FALSE
  )
  
  # top_features_v <- vector("list",nrow(tree_tb))
  # top_feature_anti_v <- vector("list",nrow(tree_tb))
  # for (i in 1:length(calc_nodes)) {
  #   ln <- calc_nodes[i]
  #   top_features_v[[ln]] <- top_features[[i]]
  #   top_feature_anti_v[[ln]] <- top_feature_anti[[i]]
  # }
  
  return(list(
    top_features=top_features, 
    top_feature_anti=top_feature_anti, 
    top_feature_sister=top_feature_sister,
    splits_mat=m_splits, heatmap=hm
  ))
}
test <- tree_gap_genes(
  tree = cooc$tree,
  feature_matrix = log2(merge_ab$ab_o + 1)
)
