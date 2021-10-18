require(futile.logger)
require(parallelDist)
require(GenomicFeatures)
require(GenomicRanges)
require(methods)
require(tidyverse)
require(reshape2)
require(cowplot)
        
infercnv_allele <- methods::setClass("infercnv_allele",
                                     slots = c(allele.data = "ANY",
                                               coverage.data = "ANY",
                                       
                                               allele.bulk.data = "ANY",
                                               coverage.bulk.data = "ANY",
                                               
                                               ref_annotation = "ANY",
                                               reference_grouped_cell_indices = "list",
                                               observation_grouped_cell_indices = "list",
                                               
                                               SNP_info = "GRanges",
                                       
                                               allele.lesser.data = "ANY",
                                               allele.lesser.bulk.data = "ANY",

                                               HMM_info = "list",
                                               HMM_region = "GRanges"))

CreateInfercnvAlleleObject <- function(allele_counts_matrix,
                                       coverage_counts_matrix,
                                       ref_annotation_file,
                                       allele_counts_bulk = NULL,
                                       coverage_counts_bulk = NULL,
                                       sample_annotation_file = NULL,
                                       ref_group_names = NULL,
                                       rowname_by = ":",
                                       delim = "\t",
                                       chr_exclude = c('chrX', 'chrY', 'chrM')){
  
  ## Read alt and tot matrix
  read_count_matrix <- function(matrix, name){
    #browser()
    if (Reduce("|", is(matrix) == "character")) {
      flog.info(sprintf("Parsing matrix: %s ...", matrix)) 
      if (substr(matrix, nchar(matrix)-2, nchar(matrix)) == ".gz") {
        raw.data <- read.table(connection <- gzfile(matrix, 'rt'), sep=delim, header=TRUE, row.names=1, check.names=FALSE) %>% as.matrix()
        close(connection)
        return(raw.data)
      }
      else if(substr(matrix, nchar(matrix)-3, nchar(matrix)) == ".rds") {
        return(readRDS(matrix))
      }
      else {
        raw.data <- read.table(matrix, sep=delim, header=TRUE, row.names=1, check.names=FALSE) %>% as.matrix()
        return(raw.data)
      }
    } 
    else if (Reduce("|", is(matrix) %in% c("dgCMatrix", "matrix"))) {
      return(matrix)
    } 
    else if (Reduce("|", is(matrix) %in% c("data.frame"))) {
      return(as.matrix(matrix))
    } 
    else {
      stop(sprintf("CreateInfercnvAlleleObject:: Error, %s isn't recognized as a matrix, data.frame, or filename",
                   name))
    }
  }
  
  allele_counts_matrix <- read_count_matrix(allele_counts_matrix, name = "allele data")
  coverage_counts_matrix <- read_count_matrix(coverage_counts_matrix, name = "coverage data")
  
  if(mean(allele_counts_matrix > 0 & allele_counts_matrix < 1) > 0){
    allele_counts_matrix[allele_counts_matrix > 0 & allele_counts_matrix < 1] <- 0 
  }
  
  if(!isTRUE(all.equal(rownames(allele_counts_matrix), rownames(coverage_counts_matrix))) |
     !isTRUE(all.equal(colnames(allele_counts_matrix), colnames(coverage_counts_matrix)))){
    stop("CreateInfercnvAlleleObject:: Error, the dimension of allele data is not same as the coverage data")
  }
  
  ## Initialize the bulk data if it's unavailable
  if(is.null(allele_counts_bulk) | is.null(coverage_counts_bulk)) {
    flog.info("Creating in-silico bulk ...")
    allele_counts_bulk <- rowSums(allele_counts_matrix > 0)
    coverage_counts_bulk <- rowSums(coverage_counts_matrix > 0)
  }
  
  ## Initialize the genome annotation
  if (Reduce("|", is(ref_annotation_file) == "character")) {
    flog.info(sprintf("Parsing gene annotation file: %s ...", ref_annotation_file)) 
    ref_annot <- read.table(ref_annotation_file, header=F, sep = delim, row.names=NULL, stringsAsFactors=F)
  }
  else if (Reduce("|", is(ref_annotation_file) %in% c("TxDb", "dgCMatrix", "matrix", "data.frame"))) {
    ref_annot <- ref_annotation_file
  }
  else {
    stop("CreateInfercnvAlleleObject:: Error, ref_annotation_file isn't recognized as a TxDb, matrix, data.frame, or filename")
  }
  if (class(ref_annot) != "TxDb"){
    colnames(ref_annot) = c('genename', 'seqnames', 'start', 'end')
  }
  
  ## Initialize cell type info
  if(!is.null(sample_annotation_file) | !is.null(ref_group_names)){
    if (Reduce("|", is(sample_annotation_file) == "character")) {
      flog.info(sprintf("Parsing cell annotation file: %s ...", sample_annotation_file))
      sample_annotation_file <- read.table(sample_annotation_file, header=FALSE, 
                                           row.names = NULL, sep = delim, stringsAsFactors = FALSE, 
                                           colClasses = c('character', 'character'))
    }
    else if (Reduce("|", is(sample_annotation_file) %in% c("dgCMatrix", "matrix", "data.frame"))) {
      sample_annotation_file <- sample_annotation_file
    }
    else {
      stop("CreateInfercnvAlleleObject:: Error, sample_annotation_file isn't recognized as a matrix, data.frame, or filename")
    }
    colnames(sample_annotation_file) <- c("cell", "celltype")
    if(mean(colnames(allele_counts_matrix) %in% sample_annotation_file$cell) != 1 |
       mean(colnames(coverage_counts_matrix) %in% sample_annotation_file$cell) != 1){
      stop("CreateInfercnvAlleleObject:: Error, some cells are missing in your cell annotation")
    }
    else{
      sample_annotation_file <- sample_annotation_file[sample_annotation_file$cell %in% colnames(allele_counts_matrix),]
      rownames(sample_annotation_file) <- sample_annotation_file$cell
      sample_annotation_file <- sample_annotation_file[colnames(allele_counts_matrix),]
    }
    
    ref_group_idx <- lapply(ref_group_names, 
                            function(x) {
                              if(sum(sample_annotation_file$celltype %in% x) > 0){
                                return(which(sample_annotation_file$celltype %in% x))
                              }
                              else
                                {stop(sprintf("Not identifying cells with classification %s", x))
                                }
                              })
    names(ref_group_idx) <- ref_group_names
    
    obs_group_names <- setdiff(sample_annotation_file$celltype, ref_group_names)
    obs_group_idx <- lapply(obs_group_names, 
                            function(x) which(sample_annotation_file$celltype %in% x))
    names(obs_group_idx) <- obs_group_names
    
  }
  else {
    ref_group_idx <- list()
    obs_group_idx <- list("Unknown" = seq_len(ncol(allele_counts_matrix)))
  }
  
  ## Initialize SNP site info
  snps.df <- rownames(allele_counts_matrix) 
  snps.df <- data.frame(do.call(rbind,strsplit(snps.df,rowname_by)), stringsAsFactors=F)
  
  if(ncol(snps.df)==2) {
    snps.df <- cbind(snps.df, snps.df[,2])
  }
  colnames(snps.df) <- c('chr','start','end')
  
  if(!grepl('chr', snps.df[1,1])) {
    snps.df[,1] <- paste0('chr', snps.df[,1])
  }
  
  ## create a Granges object for SNP
  snps <- with(snps.df, 
               GenomicRanges::GRanges(chr, 
                                      IRanges::IRanges(as.numeric(as.character(start)), 
                                                       as.numeric(as.character(end)))))
  
  names(snps) <- rownames(allele_counts_matrix) <- 
    names(allele_counts_bulk) <- 
    rownames(coverage_counts_matrix)  <- 
    names(coverage_counts_bulk) <- 
    apply(snps.df, 1, paste0, collapse=":")
  
  snps <- snps %>% sortSeqlevels() %>% sort()
  
  ## Initialize a object 
  object <- new(Class = "infercnv_allele",
                allele.data = allele_counts_matrix[names(snps),],
                coverage.data = coverage_counts_matrix[names(snps),],
                allele.bulk.data = allele_counts_bulk[names(snps)],
                coverage.bulk.data = coverage_counts_bulk[names(snps)],
                ref_annotation = ref_annot,
                reference_grouped_cell_indices = ref_group_idx,
                observation_grouped_cell_indices = obs_group_idx,
                SNP_info = snps)
  
  # validate here
  validate_infercnvallele_obj(object)
  return(object)
  
}

validate_infercnvallele_obj <- function(infercnv_allele_obj) {
  
  flog.info("Validating infercnv_allele obejct ...")
  
  if (isTRUE(all.equal(rownames(infercnv_allele_obj@allele.data), 
                       rownames(infercnv_allele_obj@coverage.data),
                       names(infercnv_allele_obj@allele.bulk.data),
                       names(infercnv_allele_obj@coverage.bulk.data),
                       names(infercnv_allele_obj@SNP_info)))) {
    return()
  }
  stop("Please check your input files carefully")
}

setAlleleMatrix <- function(infercnv_allele_obj,
                            filter = TRUE, het.deviance.threshold = 0.05, min.cell=3){
  
  if (filter){
    
    E <- infercnv_allele_obj@allele.bulk.data/infercnv_allele_obj@coverage.bulk.data
    filter_index <- E > het.deviance.threshold & E < 1-het.deviance.threshold
  
    if(sum(filter_index) < 0.01*length(infercnv_allele_obj@allele.bulk.data)) {
      flog.info("WARNING! CLONAL DELETION OR LOH POSSIBLE!")
    }
  
    infercnv_allele_obj@allele.data <- infercnv_allele_obj@allele.data[filter_index,]
    infercnv_allele_obj@coverage.data <- infercnv_allele_obj@coverage.data[filter_index,]
    infercnv_allele_obj@allele.bulk.data <- infercnv_allele_obj@allele.bulk.data[filter_index]
    infercnv_allele_obj@coverage.bulk.data <- infercnv_allele_obj@coverage.bulk.data[filter_index]
  
    filter_index <- rowSums(infercnv_allele_obj@coverage.data > 0) >= min.cell
    
    flog.info(sprintf("%s heterozygous SNPs identified ...", sum(filter_index)))
  
    infercnv_allele_obj@allele.data <- infercnv_allele_obj@allele.data[filter_index,]
    infercnv_allele_obj@coverage.data <- infercnv_allele_obj@coverage.data[filter_index,]
    infercnv_allele_obj@allele.bulk.data <- infercnv_allele_obj@allele.bulk.data[filter_index]
    infercnv_allele_obj@coverage.bulk.data <- infercnv_allele_obj@coverage.bulk.data[filter_index]
  }
  
  flog.info("Setting composite lesser allele count ...")
  
  E <- infercnv_allele_obj@allele.bulk.data/infercnv_allele_obj@coverage.bulk.data
  thr_index <- E > 0.5
  
  infercnv_allele_obj@allele.lesser.data <- infercnv_allele_obj@allele.data
  infercnv_allele_obj@allele.lesser.data[thr_index,] <- infercnv_allele_obj@coverage.data[thr_index,] - infercnv_allele_obj@allele.data[thr_index,]
  
  infercnv_allele_obj@allele.lesser.bulk.data <- infercnv_allele_obj@allele.bulk.data
  infercnv_allele_obj@allele.lesser.bulk.data[thr_index] <- infercnv_allele_obj@coverage.bulk.data[thr_index] - infercnv_allele_obj@allele.bulk.data[thr_index]
  
  #message("Mapping snps to genes ...")
  
  #gf <- ChIPseeker::annotatePeak(peak=snps, TxDb=infercnv_allele_obj@ref_annotation)
  #gf.df <- ChIPseeker:::as.data.frame.csAnno(gf)$geneId
  #snps$geneid <- gf.df
  #infercnv_allele_obj@SNP_info <- snps
  
  rownames(infercnv_allele_obj@allele.lesser.data) <- 
    names(infercnv_allele_obj@allele.lesser.bulk.data) <- 
    rownames(infercnv_allele_obj@allele.data)

  flog.info("Done setting initial allele matrices!")
  
  return(infercnv_allele_obj)
  
}

calAlleleBoundaries <- function(infercnv_allele_obj,
                                distance_method = c("Filter_threshold", "Remove_NA", "HB"), ncores = 20,
                                min.traverse = 3, t = 1e-6, pd = 0.1, pn = 0.45, 
                                min.num.snps = 5, trim = 0.1){
  
  flog.info("Starting cluster cells from population ...")
  mat.tot <- infercnv_allele_obj@allele.lesser.data/infercnv_allele_obj@coverage.data
  
  fraction_method <- match.arg(distance_method)
  flog.info(sprintf("Using %s to process matrix ...", fraction_method))
  
  if(distance_method == "Filter_threshold"){
    mat.tot[is.na(mat.tot)] <- 0 # omit no coverage
    mat.tot[infercnv_allele_obj@allele.data == 0 & infercnv_allele_obj@coverage.data != 0] <- 0.001 # pseudo count for total coverage not 0
    mat.tot[infercnv_allele_obj@coverage.data <= 2] <- 0 # filter with the cutoff of 3 coverage
    mat.smooth <- apply(mat.tot, 2, caTools::runmean, k=31)
    
    flog.info("Starting calculate distance ...")
    d <- parDist(t(mat.smooth), method = "euclidean", threads = ncores) # parallel dist
  }
  else if(distance_method == "Remove_NA"){
    mat.smooth <- apply(mat.tot, 2, caTools::runmean, k=31)
    mat.smooth[is.na(mat.smooth)] <- 0
    
    flog.info("Starting calculate distance ...")
    d <- parDist(t(mat.smooth), method = "euclidean", threads = ncores) # parallel dist
  }
  else if(distance_method == "HB"){
    mat.smooth <- apply(mat.tot, 2, caTools::runmean, k=31)
    flog.info("Starting calculate distance ...")
    d <- dist(t(mat.smooth), method = "euclidean") # too slow
    d[is.na(d)] <- 0
    d[is.nan(d)] <- 0
    d[is.infinite(d)] <- 0
  }

  flog.info("Starting calculate Hierarchical Clustering ...")
  hc <- hclust(d, method="ward.D2")
  
  flog.info('Starting iterative HMM ...')
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
    flog.info('Exiting; no new bound SNPs found ...')
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
  tbv <- tbv[tbv >= min.num.snps]
  if(length(tbv)==0) {
    flog.info(sprintf('Exiting; less than %s new bound SNPs found ...', min.num.snps))
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
  
  flog.info("Done extracting HMM regions ...")
  
  return(infercnv_allele_obj)
  
}

plot_allele <- function(infercnv_allele_obj, name_to_plot,
                        trend_smK = 31,
                        CELL_POINT_ALPHA = 0.6, dotsize=0.3, colorscheme = "BlueRed"){
  
  allele_matrix <- infercnv_allele_obj@allele.data/infercnv_allele_obj@coverage.data # fraction
  allele_matrix[is.na(allele_matrix)] <- 0 # omit no coverage
  allele_matrix[infercnv_allele_obj@allele.data == 0 & infercnv_allele_obj@coverage.data != 0] <- 0.001 # pseudo count for total coverage not 0
  allele_matrix[infercnv_allele_obj@coverage.data <= 2] <- 0 # filter with the cutoff of 3 coverage
  
  gencode_gene_pos <- infercnv_allele_obj@ref_annotation
  gencode_gene_pos$chr = str_replace(string=gencode_gene_pos$seqnames, pattern="chr", replacement="")
  gencode_gene_pos = gencode_gene_pos %>% filter(chr %in% 1:22)
  gencode_gene_pos = gencode_gene_pos %>% mutate(chr = ordered(chr, levels=1:22))
  gencode_gene_pos$pos = as.numeric( (gencode_gene_pos$start + gencode_gene_pos$end)/2)
  
  chr_maxpos = gencode_gene_pos %>% group_by(chr) %>% summarize(maxpos = max(pos))
  chr_maxpos$minpos = 1
  
  normal_cells = colnames(infercnv_allele_obj@allele.data)[unlist(infercnv_allele_obj@reference_grouped_cell_indices)]
  malignant_cells = colnames(infercnv_allele_obj@allele.data)[unlist(infercnv_allele_obj@observation_grouped_cell_indices)]
  
  num_normal_cells = length(normal_cells)
  
  flog.info("Building snp plot ...")
  
  min.cells = 3
  cells_w_ref_allele = rowSums(allele_matrix != 0 & allele_matrix < 0.5)
  cells_w_alt_allele = rowSums(allele_matrix != 0 & allele_matrix > 0.5)
  
  allele_matrix = allele_matrix[(cells_w_ref_allele >= min.cells & cells_w_alt_allele >= min.cells), ]
  
  num_snps_all = nrow(allele_matrix)
  flog.info(sprintf("Number het snps used: %s ...", num_snps_all))
  
  #malignant_cell_idx = which(colnames(allele_matrix) %in% malignant_cells)
  
  flog.info("Setting alt allele fraction to the tumor cell-population minor allele ...")
  
  mAF_allele_matrix = apply(allele_matrix, 1, function(x) {
    nonzero_val_idx = which(x>0)
    nonzero_vals = x[nonzero_val_idx]
    
    frac_high = sum(nonzero_vals>0.5)/length(nonzero_vals)
    
    ## focus allele selection based on the tumor cells only.
    tumor_vals = x[unlist(infercnv_allele_obj@observation_grouped_cell_indices)]
    tumor_nonzero_vals = tumor_vals[tumor_vals>0]
    if (length(tumor_nonzero_vals) > 0) {
      frac_high = sum(tumor_nonzero_vals>0.5)/length(tumor_nonzero_vals)
    }
    
    if ( frac_high > 0.5) {
      x[x==1] = 0.999
      x[nonzero_val_idx ] = 1 - x[nonzero_val_idx]
    }
    x
  })
  
  allele_matrix = t(mAF_allele_matrix)
  
  flog.info("Melting matrix ...")
  datamelt = melt(as.matrix(allele_matrix))
  colnames(datamelt) = c('chrpos', 'cell', 'AF')
  
  datamelt = datamelt %>% separate(chrpos, ":", into=c('seqnames', 'pos', "end"), remove=FALSE)
  
  datamelt$chr = str_replace(string=datamelt$seqnames, pattern="chr", replacement="")
  
  datamelt = datamelt %>% filter(chr %in% 1:22)
  
  datamelt = datamelt %>% mutate(chr = ordered(chr, levels=1:22))
  datamelt$pos = as.numeric(datamelt$pos)
  
  ## get chr bounds for plotting later.
  chr_maxpos = datamelt %>% group_by(chr) %>% summarize(maxpos = max(pos))
  chr_maxpos$minpos = 1
  
  
  datamelt = datamelt %>% filter(AF > 0)  ## if AF==0, then means we have no coverage.
  
  datamelt = datamelt %>% mutate(cellchr = paste(cell, seqnames, sep=":"))
  
  midpt = mean(datamelt$AF)
  
  datamelt$sample_type = "tumor"
  if (num_normal_cells > 0) {
    datamelt$sample_type[ datamelt$cell %in% normal_cells ] = "normal"
  }
  
  make_plots = function(dataToPlot, chr_maxpos, infercnv_allele_obj) {
    #browser()
    ## normal cell plot
    
    normal_snps_plot = NULL
    
    if (num_normal_cells > 0) {
      
      flog.info("Making normal plot ...")
      
      normal_dataToPlot = dataToPlot %>% filter(cell %in% normal_cells)
      
      
      normal_snps_plot = ggplot(data=normal_dataToPlot) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
        theme_bw() +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()
        ) +
        
        geom_vline(data=chr_maxpos, aes(xintercept=minpos), color=NA) +
        geom_vline(data=chr_maxpos, aes(xintercept=maxpos), color=NA) +
        
        geom_point(aes(x=pos, y=cell, color=AF), alpha=CELL_POINT_ALPHA, size=dotsize) + scale_radius()
      
      
      if (colorscheme == "BlueRed") {
        normal_snps_plot = normal_snps_plot + scale_color_gradient(low="blue", high="red")
      } else {
        normal_snps_plot = normal_snps_plot + scale_color_gradient2(low="blue", mid='yellow', high="red", midpoint=midpt)
      }
    }
    
    ## malignant cell plot
    
    flog.info("Making malignant plot ...")
    
    num_malignant_cells = length(malignant_cells)
    
    malignant_dataToPlot = dataToPlot %>% filter(cell %in% malignant_cells)
    
    malignant_snps_plot = ggplot(data=malignant_dataToPlot) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
      theme_bw() +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
      ) +
      
      geom_vline(data=chr_maxpos, aes(xintercept=minpos), color=NA) +
      geom_vline(data=chr_maxpos, aes(xintercept=maxpos), color=NA) +
      
      geom_point(aes(x=pos, y=cell, color=AF), alpha=CELL_POINT_ALPHA, size=dotsize) + scale_radius()
    
    if (colorscheme == "BlueRed") {
      malignant_snps_plot = malignant_snps_plot + scale_color_gradient(low="blue", high="red")
    } else {
      malignant_snps_plot = malignant_snps_plot + scale_color_gradient2(low="blue", mid='yellow', high="red", midpoint=midpt)
    }
    
    
    flog.info("Making allele freq plot ...")
    
    
    allele_freq_means = dataToPlot %>%
      group_by(chrpos,sample_type) %>%
      mutate(grp_pos_mean_AF = mean(AF)) %>% select(chrpos, chr, pos, sample_type, grp_pos_mean_AF) %>% unique()
    
    
    
    ## smooth mean AFs across chromosomes for each sample type.
    flog.info("Smoothing sample means for trend lines ...")
    
    allele_freq_means = allele_freq_means %>% mutate(sampleTypeChr = paste(sample_type, chr, sep=":"))
    splitdata = split(allele_freq_means, allele_freq_means$sampleTypeChr)
    
    smoother = function(df) {
      df = df %>% arrange(pos)
      df$grp_pos_mean_AF_sm = caTools::runmean(df$grp_pos_mean_AF, k=trend_smK, align="center")
      return(df)
    }
    
    allele_freq_means = do.call(rbind, lapply(splitdata, smoother))
    
    allele_freq_plot = allele_freq_means %>%
      ggplot(aes(x=pos, y=grp_pos_mean_AF)) +
      facet_grid (~chr, scales = 'free_x', space = 'fixed') +
      theme_bw() +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
      ) +
      geom_vline(data=chr_maxpos, aes(xintercept=minpos), color=NA) +
      geom_vline(data=chr_maxpos, aes(xintercept=maxpos), color=NA) +
      geom_point(aes(color=sample_type), alpha=0.2, size=0.2)
    
    
    allele_freq_plot_w_trendlines = allele_freq_plot +
      geom_line(data=allele_freq_means,
                aes(x=pos, y=grp_pos_mean_AF_sm, color=sample_type), size=0.5, alpha=1)
    
    # HMM prediction if existed
    if(!is.null(infercnv_allele_obj@HMM_region)){
      
      flog.info("Making HMM prediction plot ...")
      #HMM_candidate <- infercnv_allele_obj@HMM_region[seqnames(infercnv_allele_obj@HMM_region) %in% paste0("chr", 1:22),]
      malignant_dataToPlot_Gr <- with(malignant_dataToPlot, 
                                      GenomicRanges::GRanges(seqnames, 
                                                             IRanges::IRanges(as.numeric(as.character(pos)), 
                                                                              as.numeric(as.character(end)))))
      malignant_dataToPlot$HMM_region <- "Neural"
      malignant_dataToPlot$HMM_region[malignant_dataToPlot_Gr %over% infercnv_allele_obj@HMM_region] <- "Candidate Deletions/LOHs"
      
      malignant_HMM_plot = ggplot(data=malignant_dataToPlot) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
        theme_bw() +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()
        ) +
        geom_vline(data=chr_maxpos, aes(xintercept=minpos), color=NA) +
        geom_vline(data=chr_maxpos, aes(xintercept=maxpos), color=NA) +
        geom_point(aes(x=pos, y=cell, color=HMM_region), alpha=CELL_POINT_ALPHA, size=dotsize) + scale_radius()
    }
    
    if (num_normal_cells > 0) {
      
      ratio_normal_cells = max(0.25, num_normal_cells/(num_normal_cells + num_malignant_cells))
      
      pg = plot_grid(normal_snps_plot, allele_freq_plot_w_trendlines, 
                     malignant_HMM_plot, malignant_snps_plot, 
                     ncol=1, align='v', rel_heights=c(0.25, 0.25, 0.25, 0.25))
      
      return(pg)
      
    } else {
      pg = plot_grid(allele_freq_plot_w_trendlines, malignant_HMM_plot, malignant_snps_plot, 
                     ncol=1, align='v', rel_heights=c(0.25, 0.25, 0.5))
      
      return(pg)
    }
    
  }
  
  flog.info("Generating outputs ...")
  
  p <- make_plots(datamelt, chr_maxpos, infercnv_allele_obj)
  ggsave (name_to_plot, p, width = 13.33, height = 7.5, units = 'in', dpi = 300)
  flog.info("Done!")

}