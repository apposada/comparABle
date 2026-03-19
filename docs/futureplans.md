# Future plans for comparABle

* Make the code more modular: smaller functions for specific tasks
  * This includes making the plot functions less clunky and more modular too.
* Make the code smarter: detect what parameters and objects have been passed to silently run (or not) certain analyses (able to handle what to run if some attribute files are not specified). I.e. if no GO data is passed down, do not run GO analyses.
* Consider constructing "objects" where to store the data, that the functions can detect.
* Add fast, robust, and reproducible functions for:
  * Normalisation (quantile normalisation, scaling, log transformation);
  * Correlation (consider `WGCNA::cor()` or `tgstat::cor()`, for example), Jensen-Shannon Divergence;
  * Additional metrics for divergence/similarity: Jaccard similarity between pairs of modules
  * Regression/fitting for data with heteroscedasticity (see [this thread](https://support.bioconductor.org/p/9151444/));
* Add additional methods:
  * Consider the case of treeExp (see [original paper](https://doi.org/10.1093/bioinformatics/btz405), [Github](https://github.com/hr1912/TreeExp), [tutorial](https://jingwyang.github.io/TreeExp-Tutorial/) );
  * Iterative Comparison of Coexpression ([original paper](https://link.springer.com/article/10.1186/gb-2007-8-4-r50), [implementation by Sebé-Pedrós lab](https://github.com/sebepedroslab/metacell-downstream-functions/blob/master/Cross_species_functions.R#L3196))
  * Additional ways to identify genes. Currently it is limma-based, this works but there might be other ways. Maybe it could be later expanded into an aggregated score.
  * Related to this last bit: add the second part of `treefromEnsemblCluster()` found in [the original repo](ttps://github.com/sebepedroslab/Stylophora_single_cell_atlas/blob/857eb758bb6886bd91482cfe601e9bd5f56b12de/metacell_downstream_functions/Tree_functions.R);
  * For WGCNA, multi-species: method to identify stable sub-circuits of co-enriched modules with a species composition representative of the phylogeny at hand (might require defining a tree);
  * 
* Expand highly-correlating gene heuristics to **groups** of samples and not only **pairs** of samples (using graphs (derived from JSD/similarity matrices) + communities, half-baked already in ongoing project);
* Add proper documentation;
* Future future plans: add support for multiple species simultaneously;
* Make a cuter logo. This should be higher in the list of priorities.
