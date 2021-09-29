##### This script aims to extract potential regions of LOH events based on the allele info using HMM model

## packages needed
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

## defined functions
setAlleleMatrix <- function(alter_sc, cov_sc,
                            alter_bulk = NULL, cov_bulk = NULL,
                            filter = TRUE, het.deviance.threshold = 0.05, min.cell=3,
                            verbose = TRUE){
  
  if(verbose){
    cat("Initializing allele matrices ... \n")
  }
  
  if(is.null(alter_bulk) | is.null(cov_bulk)) {
    
    if(verbose) {
      cat("Creating in-silico bulk ... \n")
      cat(paste0("using ", ncol(alter_sc), " cells ... \n"))
    }
    
    alter_bulk <- rowSums(alter_sc > 0)
    cov_bulk <- rowSums(cov_sc > 0)
    
  } 
  
  if(filter) {
    
    if(verbose) {
      cat("Filtering for putative heterozygous snps ... \n")
      cat(paste0("allowing for a ", het.deviance.threshold, 
                 " deviation from the expected 0.5 heterozygous allele fraction ... \n"))
    }
    
    E <- alter_bulk/cov_bulk
    vi <- E > het.deviance.threshold & E < 1-het.deviance.threshold
    
    if(verbose) {
      cat(paste0(sum(vi), " heterozygous SNPs identified \n"))
    }
    
    if(sum(vi) < 0.01*length(alter_bulk)) {
      cat("WARNING! CLONAL DELETION OR LOH POSSIBLE! \n")
    }
    
    alter_sc <- alter_sc[vi,]
    cov_sc <- cov_sc[vi,]
    alter_bulk <- alter_bulk[vi]
    cov_bulk <- cov_bulk[vi]
    
    ## must have coverage in at least n cells
    if(verbose) {
      cat(paste0("must have coverage in at least ", min.cell, " cells ... \n"))
    }
    
    vi <- rowSums(cov_sc > 0) >= min.cell
    cat(paste0(sum(vi), " heterozygous SNPs identified \n"))
    alter_sc <- alter_sc[vi,]
    cov_sc <- cov_sc[vi,]
    alter_bulk <- alter_bulk[vi]
    cov_bulk <- cov_bulk[vi]
    
  }
  
  if(verbose) {
    cat("Setting composite lesser allele count ... \n")
  }
  
  vi <- alter_bulk/cov_bulk > 0.5
  alter_less_sc <- alter_sc
  alter_less_sc[vi,] <- cov_sc[vi,] - alter_sc[vi,]
  
  alter_less_bulk <- alter_bulk
  alter_less_bulk[vi] <- cov_bulk[vi] - alter_bulk[vi]
  
  if(verbose) {
    cat("Setting SNP regions ... \n")
  }
  
  snps.df <- rownames(alter_sc)
  snps.df <- data.frame(do.call(rbind,strsplit(snps.df,":|-| ")), stringsAsFactors=F)
  
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
  names(snps) <- rownames(alter_less_sc) <- names(alter_less_bulk) <- rownames(cov_sc) <- names(cov_bulk) <- apply(snps.df, 1, paste0, collapse=":")
  
  if(verbose) {
    cat("Done setting initial allele matrices! \n")
  }
  
  return(list("alter_less_sc" = alter_less_sc,
              "alter_less_bulk" = alter_less_bulk,
              "cov_sc" = cov_sc,
              "cov_bulk" = cov_bulk,
              "snps" = snps))
  
}

setGeneMapping <- function(site, txdb,
                           fill=TRUE, verbose=TRUE){
  
  if(verbose) {
    cat("Mapping snps to genes ... \n")
  }
  
  gf <- ChIPseeker::annotatePeak(peak=site, TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
  gf.df <- ChIPseeker:::as.data.frame.csAnno(gf)$geneId
  names(gf.df) <- names(site)
  
  if(!fill) {
    gf.df[is.na(ChIPseeker:::as.data.frame.csAnno(gf)$annotation)] <- NA
    gf.df <- na.omit(gf.df)
  }
  
  geneFactor <- gf.df
  if(verbose) {
    cat("Done mapping snps to genes! \n")
  }
  
  return(geneFactor)

}

calAlleleBoundaries <- function(r.sub, n.sc.sub,
                                min.traverse = 3, t = 1e-6, pd = 0.1, pn = 0.45, 
                                min.num.snps = 5, trim = 0.1, verbose = TRUE){
  
  pred.snps.r <- matrix(0, nrow(r.sub), ncol(r.sub))
  rownames(pred.snps.r) <- rownames(r.sub)
  colnames(pred.snps.r) <- colnames(r.sub)
  bound.snps.final <- list()
  
  ## lesser allele fraction
  mat.tot <- r.sub/n.sc.sub
  
  mat.smooth <- apply(mat.tot, 2, caTools::runmean, k=31)
  d <- dist(t(mat.smooth))
  d[is.na(d)] <- 0
  d[is.infinite(d)] <- 0
  hc <- hclust(d, method="ward.D2")
  
  if(verbose) {
    cat('iterative HMM ... ')
  }
  
  ## iterative HMM
  heights <- 1:min(min.traverse, ncol(r.sub))
  ## cut tree at various heights to establish groups
  boundsnps.pred <- lapply(heights, function(h) {
    
    ct <- cutree(hc, k = h)
    cuts <- unique(ct)
    
    ## look at each group, if deletion present
    boundsnps.pred <- lapply(cuts, function(group) {
      
      if(sum(ct == group)>1) {
        mafl <- rowSums(r.sub[, ct == group]>0)
        sizel <- rowSums(n.sc.sub[, ct == group]>0)
        
        ## change point
        delta <- c(0, 1)
        z <- HiddenMarkov::dthmm(mafl, matrix(c(1-t, t, t, 1-t), 
                                              byrow=TRUE, nrow=2), 
                                 delta, "binom", list(prob=c(pd, pn)), 
                                 list(size=sizel), discrete=TRUE)
        results <- HiddenMarkov::Viterbi(z)
        
        ## Get boundaries from states
        boundsnps <- rownames(r.sub)[results == 1]
        return(boundsnps)
        
      }
    })
  })
  
  boundsnps_res <- table(unlist(boundsnps.pred))
  
  ## vote
  vote <- rep(0, nrow(r.sub))
  names(vote) <- rownames(r.sub)
  vote[names(boundsnps_res)] <- boundsnps_res
  
  if(verbose) {
    cat(paste0('max vote:', max(vote), '\n\n'))
  }
  
  if(max(vote) == 0) {
    if(verbose) {
      cat('Exiting; no new bound SNPs found.\n')
    }
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
  
  HMM_region <- lapply(names(tbv), function(ti) {
    
    bound.snps.new <- names(bound.snps.cont)[bound.snps.cont == ti]
    
    ## trim
    bound.snps.new <- bound.snps.new[1:round(length(bound.snps.new)-length(bound.snps.new)*trim)]
    
    if(verbose) {
      cat('SNPS AFFECTED BY DELETION/LOH: \n')
      cat(bound.snps.new)
      cat("\n\n")
    }
    
    return(bound.snps.new)
    
  })
    
  return(HMM_region)

}

plotAlleleProfile <- function(r.sub, n.sc.sub, 
                              l.sub, n.bulk.sub, snps,
                              region=NULL, chrs=paste0('chr', c(1:22)), 
                              setWidths=FALSE, cellOrder=NULL, 
                              filter=FALSE, returnPlot=FALSE, max.ps=3){
  
  if(!is.null(region)) {
    overlap <- IRanges::findOverlaps(region, snps)
    ## which of the ranges did the position hit
    hit <- rep(FALSE, length(snps))
    hit[S4Vectors::subjectHits(overlap)] <- TRUE
    if(sum(hit) < 10) {
      cat(paste0("WARNING! ONLY ", sum(hit), " SNPS IN REGION! \n"))
    }
    vi <- hit
    r.sub <- r.sub[vi,]
    n.sc.sub <- n.sc.sub[vi,]
    l.sub <- l.sub[vi]
    n.bulk.sub <- n.bulk.sub[vi]
    
    chrs <- region@seqnames@values
  }
  if(!is.null(cellOrder)) {
    r.sub <- r.sub[,cellOrder]
    n.sc.sub <- n.sc.sub[,cellOrder]
  }
  if(filter) {
    ## filter out snps without coverage
    vi <- rowSums(n.sc.sub) > 0
    r.sub <- r.sub[vi,]
    n.sc.sub <- n.sc.sub[vi,]
    l.sub <- l.sub[vi]
    n.bulk.sub <- n.bulk.sub[vi]
  }
  
  r.tot <- cbind(r.sub/n.sc.sub, 'Bulk'=l.sub/n.bulk.sub)
  n.tot <- cbind(n.sc.sub, 'Bulk'=n.bulk.sub)
  
  require(ggplot2)
  require(reshape2)
  
  if(setWidths) {
    widths <- sapply(chrs, function(chr) {
      sum(grepl(paste0('^',chr,':'), rownames(r)))
    }); widths <- widths/max(widths)*100
  } else {
    widths <- rep(1, length(chrs))
  }
  
  plist <- lapply(chrs, function(chr) {
    vi <- grepl(paste0('^',chr,':'), rownames(r.tot))
    m <- melt(t(r.tot[vi,]))
    colnames(m) <- c('cell', 'snp', 'alt.frac')
    rownames(m) <- paste(m$cell, m$snp)
    m$alt.frac[is.nan(m$alt.frac)] <- NA
    n <- melt(t(n.tot[vi,]))
    colnames(n) <- c('cell', 'snp', 'coverage')
    rownames(n) <- paste(n$cell, n$snp)
    n$coverage[n$coverage>30] <- 30  # max for visualization purposes
    ##n$coverage <- log10(n$coverage+1)
    n$coverage <- n$coverage^(1/3) # cube root for visualization purposes only
    n$coverage[n$coverage==0] <- NA # if no coverage, just don't show
    dat <- cbind(m, coverage=n$coverage)
    
    p <- ggplot(dat, aes(snp, cell)) +
      ## geom_tile(alpha=0) +
      geom_point(aes(colour = alt.frac, size = coverage), na.rm=TRUE) +
      scale_size_continuous(range = c(0, max.ps)) +
      ## scale_colour_gradientn(colours = rainbow(10)) +
      scale_colour_gradient2(mid="yellow", low = "turquoise", high = "red", midpoint=0.5) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none",
        plot.margin=unit(c(0,0,0,0), "cm"),
        panel.border = element_rect(fill = NA, linetype = "solid", colour = "black"),
        plot.title = element_text(hjust = 0.5)
      ) + labs(title = chr)
    ## theme(
    ##     ## axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=rel(0.5),lineheight=1),
    ##     ## axis.text.y=element_blank(),
    ##     axis.title.y=element_blank(),
    ##     axis.ticks.y=element_blank(),
    ##     ##axis.text.y=element_text(size=rel(0.5))
    ##     legend.position="bottom"
    ##     ##panel.margin=unit(0 , "lines")
    ## )
    return(p)
  })
  
  if(returnPlot) {
    return(plist)
  } else {
    require(gridExtra)
    do.call("grid.arrange", c(plist, ncol=length(plist)))
    ##print(p)
  }

}