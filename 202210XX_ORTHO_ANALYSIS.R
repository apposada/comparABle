require(data.table)
require(stringr)
require(tidytree)
require(ggtree)
require(rvcheck)
require(treeio)
require(dendextend)
require(phangorn)
require(phytools)
require(ComplexHeatmap)
require(foreach)
require(doParallel)
require(colorspace)

a = pfla_emblar_allexpr
b = unique(RNA_blan_emblar); rownames(b) <- b$V1
o = unique(oma_blan_pfla)
f = intid_omaid[grep("TCONS|BRALA",intid_omaid$id),]
f$id = gsub("\\.p.*","",f$id)
ma = data.frame(
  id = RNA_pfla_emblar$id,
  module = RNA_pfla_emblar$newcl
)
mb = merge(
  RNA_blan_emblar[
    ,
    c("Row.names","V1")
  ]
  ,
  data.frame(
    id=rownames(RNA_blan),
    module=RNA_blan$newcl
  ),
  by.x=1,
  by.y=1
)[,c(2,3)]


samples_a = levels(condition_x)
a = a[,!(sapply(a, is.character))]
a = qnorm(a)
a = tidyup(a, highlyvariable = T)
a = rep2means(samples_a,a)

b = b[,!(sapply(b, is.character))]
samples_b = unique(sub("_.$", "", colnames(b)))
b = qnorm(b)
b = tidyup(b, highlyvariable = T)
b = rep2means(samples_b,b)

o = pair2id(o)


colnames(ma) <- c("id","module")
colnames(mb) <- c("id","module")

merge_ab <- mergedata(a,b,o)
cors <- rawcorsp(merge_ab$a_o,merge_ab$b_o)
pi <- prcomp(t(merge_ab$ab_o))

set.seed(4343)
h <- c(0.75,0.95)
clustering_algorithm <- "hclust"
clustering_method <- "average"
cor_method <- "pearson"
p <- 0.1
vargenes = rownames(merge_ab$ab_o)

cooc <- treeFromEnsembleClustering(
  x=merge_ab$ab_o, p=p, h=h,  n = 1000, vargenes = rownames(merge_ab$ab_o), bootstrap=FALSE,
  clustering_algorithm=clustering_algorithm, clustering_method=clustering_method, 
  cor_method=cor_method
)

# Co-occurrence matrix
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


plot_pca_ab(pca = pi,ab_o = merge_ab)
Heatmap(
  cors$sp,
  name = "correlation",
  cluster_rows = F,
  cluster_columns = F,
  col = rev(sequential_hcl(10,"YlOrRd"))
  )
heatmap(
  cooc$cooccurrence,
  symm = T,
  Rowv = cooc$tree$edge
  )


# function modules
modulecomp_ab <- comparemodules(ma,mb,f)

pheatmap::pheatmap(
  modulecomp_ab$logbinom,
  col = rev(sequential_hcl(10,"Blues 3")),
  cluster_cols = F,
  cluster_rows = F
  )

pheatmap::pheatmap(
  modulecomp_ab$loghypg,
  col = rev(sequential_hcl(10,"YlOrBr")),
  cluster_cols = F,
  cluster_rows = F
  )

test <- commongenes_cor(
  merge_ab$ab_o,
  across_a = c("05_MG","06_MG","07_LG","08_To"),
  across_b = c("10_27h","11_36h","12_50h","13_60h"),
  min_cor = 0.9
  )


