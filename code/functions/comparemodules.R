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
comparemodules <- function(ma, mb, f, cutoff = 30, verbose = FALSE){
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
  
  #' Create lists of gfams per module
  message('extracting families for modules from a')
  ma_f <- lapply(ma_list,lookup, d = f)
  message('extracting families for modules from b')
  mb_f <- lapply(mb_list,lookup, d = f)
  
  #' Create a matrix to store results (-logpval or similar)
  PVS <- data.frame() # create matrix to store result pval
  PVS_binom <- data.frame()
  
  
  #' Compute the upper tail of the hypergeometric 
  #' distribution (survival function) for each pair of modules
  #'  as a metric of enrichment #comment from panos @ skarmetalab
  for(i in names(ma_f)){
    a_fams <- ma_f[[i]]
    i_ <- which(names(ma_f) == i)
    if(verbose) message(i)
    
    for(j in names(mb_f)){
      b_fams <- mb_f[[j]]
      j_ <- which(names(mb_f) == j)
      if(verbose) message(j)
      
      # HYPERGEOMETRIC TEST
      sample_size <- length(a_fams)
      success_in_pop <- length(b_fams)
      common_fams <- paste(
        a_fams[which(a_fams %in% b_fams)], collapse = ","
      )
      exclusive_fams_a <- paste(
        a_fams[which(!(a_fams %in% b_fams))], collapse = ","
      )
      exclusive_fams_b <- paste(
        b_fams[which(!(b_fams %in% a_fams))], collapse = ","
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
            i, j, 
            success_in_sample, sample_size, 
            success_in_pop, gene_POP, as.numeric(hypg), 
            -log(hypg), as.numeric(binom), common_fams, # add here which are the names of the gfams enriched.
            exclusive_fams_a, # exclusive fams in a
            exclusive_fams_b # exclusive fams in b
          )
        )
      } else {
        hypg <- 1
        binom <- 1
      }
      
      PVS[i_, j_] <- hypg
      PVS_binom[i_, j_] <- binom
      
    }
  }
  
  # Tidy up of data 
  rownames(PVS) <- names(ma_list)
  colnames(PVS) <- names(mb_list)
  rownames(PVS_binom) <- names(ma_list)
  colnames(PVS_binom) <- names(mb_list)
  PVS <- as.matrix(PVS)
  PVS_binom <- as.matrix(PVS_binom)
  
  fon <- fon[!(fon$module_a == "none"), ]
  fon$hypgeom_pval <- as.numeric(fon$hypgeom_pval)
  fon$hypgeom_log <- as.numeric(fon$hypgeom_log)
  fon$binom_pval <- as.numeric(fon$binom_pval)
  
  loghypg <- -log(PVS)
  loghypg[loghypg > cutoff] <- cutoff
  
  logbinom <- -log(PVS_binom)
  logbinom[logbinom > cutoff] <- cutoff
  
  
  res <- list(
    stats = fon, 
    hypgeom = PVS, 
    binom = PVS_binom, 
    loghypg = loghypg, 
    logbinom = logbinom
  )
  
  return(res)
}
