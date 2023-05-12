# comparABle🦄🔁🐉

![logo](/graphics/comparable_logo.png?raw=true)

comparABle is a small collection of functions in R to **analyse and compare gene expression data across pairs of species**.


# Requirements

* `a`: a dataset of genes x condition for species A
* `b`: a dataset of genes x condition for species B
* `o`: a 1-to-1 orthologue association file between genes of species A and B (e.g. blast reciprocal hits, oma 1:1 orthologues, ...)
* `f`: a gene family file associating genes from species A and B (e.g. Orthofinder output
* (optional, reccommended) `ma, mb`: gene classification system. For each species a and b, a file associating genes to gene 
moules (e.g. k-means clusters, hierarchical clusters, WGCNA modules, etc.)
* (optional, reccommended) `g`: a gene age file for the gene families, and the genes, of each species

Briefly, comparABle generates correlation matrices across conditions of species, retrieves the most similar
conditions based on bootstrapped survival clustering, and compares gene modules through orthology overlap strategy:

 1. Normalise the expression datasets
 2. Associate the datasets using the association file of 1-to-1 orthologues
 3. Perform a principal component analysis with the merge of the two datasets
 4. Perform pairwise correlations and distance metrics between conditions (in our case, developmental stages) of datasets to see their relationships of similarity. By default: Pearson Correlation, Spearman Correlation, and Jensen-Shannon Divergence. 
 5. Identify the genes responsible for the pairwise-similarities observed in step4 across conditions (in our case, stages), and perform Gene Age and Gene Ontology enrichment analysis
 6. Perform a co-occurrence analysis in the merge of the two datasets (pairwise correlations, hierarchical clustering, and bootstrapping with downsampling) (Current implementation as seen in Levy et al., 2021 repo: https://github.com/sebepedroslab/Stylophora_single_cell_atlas/blob/857eb758bb6886bd91482cfe601e9bd5f56b12de/metacell_downstream_functions/Tree_functions.R )
 7. Perform Orthology Overlap analysis (hypergeometric and binomial) in pairwise correlations between gene sets (==stage-specific clusters) of each species.
 8. Identify the genes of the gene families that appear significantly enriched across pairs of gene sets, and perform Gene Age, Gene Ontology, and Functional Category enrichment analyses.


Optional analyses include downstream GO, COG functional category, and Gene Age enrichment of common genes across datasets, or common gene families across modules of the species.

comparABle can be passed the **same species in a and b** to generate an **intra-specific comparison** of conditions, and gene modules. Potentially useful when denovo-exploring novel cell types or structures at the molecular level.

In the future, I aim for comparABle to be more modular (implemented in smaller pieces of code that can be run independently) and a bit smarter (able to handle what to run if some attribute files are not specified), and more metrics will be added, as well as better ways to identify key genes.
