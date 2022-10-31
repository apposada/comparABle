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
conditions based on bootstrapped survival clustering, and compares gene modules through orthology overlap strategy.

Optional analyses include downstream GO, COG functional category, and Gene Age enrichment of common genes across datasets, or common gene families across modules of the species.

comparABle can be passed the **same species in a and b** to generate an **intra-specific comparison** of conditions, and gene modules. Potentially useful when denovo-exploring novel cell types or structures at the molecular level.
