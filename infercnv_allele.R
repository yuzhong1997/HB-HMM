infercnv_allele <- methods::setClass("infercnv_allele",
                                     slots = c(allele.lesser.data = "ANY",
                                               allele.lesser.bulk.data = "ANY",
                                               
                                               SNP_info = "ANY",
                                               
                                               allele.data = "ANY",
                                               coverage.data = "ANY",
                                               
                                               allele.bulk.data = "ANY",
                                               coverage.bulk.data = "ANY",
                                               
                                               ref_annotation = "ANY",
                                               reference_grouped_cell_indices = "list",
                                               observation_grouped_cell_indices = "list",
                                               
                                               HMM_info = "ANY",
                                               HMM_region = "ANY"))

CreateInfercnvAlleleObject <- function(allele_counts_matrix,
                                       coverage_counts_matrix,
                                       ref_annotation_file,
                                       allele_counts_bulk = NULL,
                                       coverage_counts_bulk = NULL,
                                       sample_annotation_file = NULL,
                                       ref_group_names = NULL,
                                       delim = "\t",
                                       chr_exclude = c('chrX', 'chrY', 'chrM')){
  
  if(is.null(allele_counts_bulk) | is.null(coverage_counts_bulk)) {
    
    message("No bulk datasets provided ...\nCreating in-silico bulk ...")
    
    allele_bulk <- rowSums(allele_counts_matrix > 0)
    coverage_bulk <- rowSums(coverage_counts_matrix > 0)
    
  }
  
  object <- new(Class = "infercnv_allele",
                allele.data = allele_counts_matrix,
                coverage.data = coverage_counts_matrix,
                allele.bulk.data = allele_bulk,
                coverage.bulk.data = coverage_bulk,
                ref_annotation = ref_annotation_file)
  
  return(object)
  
}

setAlleleMatrix <- function(infercnv_allele_obj,
                            filter = TRUE, het.deviance.threshold = 0.05, min.cell=3){
  
  E <- infercnv_allele_obj@allele.bulk.data/infercnv_allele_obj@coverage.bulk.data
  filter_index <- E > het.deviance.threshold & E < 1-het.deviance.threshold
  
  if(sum(filter_index) < 0.01*length(infercnv_allele_obj@allele.bulk.data)) {
    message("WARNING! CLONAL DELETION OR LOH POSSIBLE!")
  }
  
  infercnv_allele_obj@allele.data <- infercnv_allele_obj@allele.data[filter_index,]
  infercnv_allele_obj@coverage.data <- infercnv_allele_obj@coverage.data[filter_index,]
  infercnv_allele_obj@allele.bulk.data <- infercnv_allele_obj@allele.bulk.data[filter_index]
  infercnv_allele_obj@coverage.bulk.data <- infercnv_allele_obj@coverage.bulk.data[filter_index]
  
  filter_index <- rowSums(infercnv_allele_obj@coverage.data > 0) >= min.cell
  
  message(paste0(sum(filter_index), " heterozygous SNPs identified ..."))
  
  infercnv_allele_obj@allele.data <- infercnv_allele_obj@allele.data[filter_index,]
  infercnv_allele_obj@coverage.data <- infercnv_allele_obj@coverage.data[filter_index,]
  infercnv_allele_obj@allele.bulk.data <- infercnv_allele_obj@allele.bulk.data[filter_index]
  infercnv_allele_obj@coverage.bulk.data <- infercnv_allele_obj@coverage.bulk.data[filter_index]
  
  message("Setting composite lesser allele count ...")
  
  thr_index <- infercnv_allele_obj@allele.bulk.data/infercnv_allele_obj@coverage.bulk.data > 0.5
  infercnv_allele_obj@allele.lesser.data <- infercnv_allele_obj@allele.data
  infercnv_allele_obj@allele.lesser.data[thr_index,] <- infercnv_allele_obj@coverage.data[thr_index,] - infercnv_allele_obj@allele.data[thr_index,]
  
  infercnv_allele_obj@allele.lesser.bulk.data <- infercnv_allele_obj@allele.bulk.data
  infercnv_allele_obj@allele.lesser.bulk.data[thr_index] <- infercnv_allele_obj@coverage.bulk.data[thr_index] - infercnv_allele_obj@allele.bulk.data[thr_index]
  
  snps.df <- rownames(infercnv_allele_obj@allele.data)
  snps.df <- data.frame(do.call(rbind,strsplit(snps.df,":")), stringsAsFactors=F)
  
  if(ncol(snps.df)==2) {
    snps.df <- cbind(snps.df, snps.df[,2])
  }
  colnames(snps.df) <- c('chr','start','end')
  
  if(!grepl('chr', snps.df[1,1])) {
    snps.df[,1] <- paste0('chr', snps.df[,1])
  }
  
  snps <- with(snps.df, 
               GenomicRanges::GRanges(chr, 
                                      IRanges::IRanges(as.numeric(as.character(start)), 
                                                       as.numeric(as.character(end)))))
  
  message("Mapping snps to genes ...")
  
  gf <- ChIPseeker::annotatePeak(peak=snps, TxDb=infercnv_allele_obj@ref_annotation)
  gf.df <- ChIPseeker:::as.data.frame.csAnno(gf)$geneId
  snps$geneid <- gf.df
  infercnv_allele_obj@SNP_info <- snps
  
  names(infercnv_allele_obj@SNP_info) <- rownames(infercnv_allele_obj@allele.data) <- 
    names(infercnv_allele_obj@allele.bulk.data) <- 
    rownames(infercnv_allele_obj@allele.lesser.data) <- 
    names(infercnv_allele_obj@allele.lesser.bulk.data) <- 
    rownames(infercnv_allele_obj@coverage.data)  <- 
    names(infercnv_allele_obj@coverage.bulk.data) <- 
    apply(snps.df, 1, paste0, collapse=":")
  
  message("Done setting initial allele matrices!")
  
  return(infercnv_allele_obj)
  
}

calAlleleBoundaries <- function(infercnv_allele_obj,
                                min.traverse = 3, t = 1e-6, pd = 0.1, pn = 0.45, 
                                min.num.snps = 5, trim = 0.1){
  
  mat.tot <- infercnv_allele_obj@allele.lesser.data/infercnv_allele_obj@coverage.data
  
  mat.smooth <- apply(mat.tot, 2, caTools::runmean, k=31)
  d <- dist(t(mat.smooth))
  d[is.na(d)] <- 0
  d[is.infinite(d)] <- 0
  hc <- hclust(d, method="ward.D2")
  
  message('iterative HMM ... ')
  
  heights <- 1:min(min.traverse, ncol(infercnv_allele_obj@allele.lesser.data))
  
  boundsnps.pred <- lapply(heights, function(h) {
    
    ct <- cutree(hc, k = h)
    cuts <- unique(ct)
    
    ## look at each group, if deletion present
    boundsnps.pred <- lapply(cuts, function(group) {
      
      if(sum(ct == group)>1) {
        mafl <- rowSums(infercnv_allele_obj@allele.lesser.data[, ct == group]>0)
        sizel <- rowSums(infercnv_allele_obj@coverage.data[, ct == group]>0)
        
        ## change point
        delta <- c(0, 1)
        z <- HiddenMarkov::dthmm(mafl, matrix(c(1-t, t, t, 1-t), 
                                              byrow=TRUE, nrow=2), 
                                 delta, "binom", list(prob=c(pd, pn)), 
                                 list(size=sizel), discrete=TRUE)
        results <- HiddenMarkov::Viterbi(z)
        
        ## Get boundaries from states
        boundsnps <- rownames(infercnv_allele_obj@allele.lesser.data)[results == 1]
        return(boundsnps)
        
      }
    })
  })
  
  boundsnps_res <- table(unlist(boundsnps.pred))
  
  ## vote
  vote <- rep(0, nrow(infercnv_allele_obj@allele.lesser.data))
  names(vote) <- rownames(infercnv_allele_obj@allele.lesser.data)
  vote[names(boundsnps_res)] <- boundsnps_res
  
  if(max(vote) == 0) {
    
    message('Exiting; no new bound SNPs found.\n')
    
    return() ## exit iteration, no more bound SNPs found
    
  }
    
  vote[vote > 0] <- 1
  mv <- 1 ## at least 1 vote
  cs <- 1
  bound.snps.cont <- rep(0, length(vote))
  names(bound.snps.cont) <- names(vote)
  
  for(i in 2:length(vote)) {
    if(vote[i] >= mv & vote[i] == vote[i-1]) {
      bound.snps.cont[i] <- cs
    } else {
      cs <- cs + 1
    }
  }
  
  tb <- table(bound.snps.cont)
  tbv <- as.vector(tb)
  names(tbv) <- names(tb)
  tbv <- tbv[-1] # get rid of 0
  
  ## all detected deletions have fewer than 5 SNPs...reached the end
  tbv[tbv < min.num.snps] <- NA
  tbv <- na.omit(tbv)
  if(length(tbv)==0) {
    if(verbose) {
      cat(paste0('Exiting; less than ', min.num.snps, ' new bound SNPs found.\n'))
    }
    return()
  }
  
  HMM_info <- lapply(names(tbv), function(ti) {
    
    bound.snps.new <- names(bound.snps.cont)[bound.snps.cont == ti]
    
    ## trim
    bound.snps.new <- bound.snps.new[1:round(length(bound.snps.new)-length(bound.snps.new)*trim)]
    
    #if(verbose) {
    #  cat('SNPS AFFECTED BY DELETION/LOH: \n')
    #  cat(bound.snps.new)
    #  cat("\n\n")
    #}
    
    return(bound.snps.new)
    
  })
  HMM_region <- do.call("c", lapply(HMM_info, function(bs) range(infercnv_allele_obj@SNP_info[bs])))
  
  infercnv_allele_obj@HMM_info <- HMM_info
  infercnv_allele_obj@HMM_region <- HMM_region
  
  message("Done extracting HMM regions ...")
  
  return(infercnv_allele_obj)
  
}