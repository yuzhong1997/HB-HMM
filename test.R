setwd("D:/Dropbox/Broad/Script/HMM_dir")

library(magrittr)
alt.table <- read.table("sample_data/example.alt.matrix")
alt.table <- apply(alt.table, c(1,2), function(x) ifelse(x > 0 & x < 1, 0, x)) 
mean(alt.table > 0 & alt.table < 1)
dim(alt.table)

rownames(alt.table) <- gsub("::", ":", rownames(alt.table))

tot.table <- read.table("sample_data/example.tot.matrix") %>% as.matrix()
mean(tot.table > 0 & tot.table < 1)
dim(tot.table)

rownames(tot.table) <- gsub("::", ":", rownames(tot.table))

source("HMM_code.R")
allele_matrix <- setAlleleMatrix(alter_sc = alt.table, cov_sc = tot.table, 
                                 het.deviance.threshold = 0.1, verbose = TRUE)
allele_HMM <- calAlleleBoundaries(r.sub = allele_matrix$alter_less_sc, 
                                  n.sc.sub = allele_matrix$cov_sc,
                                  snps = allele_matrix$snps)
allele_HMM$Region

source("infercnv_allele.R")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
test <- CreateInfercnvAlleleObject(allele_counts_matrix = alt.table,
                                   coverage_counts_matrix = tot.table,
                                   ref_annotation_file = txdb)
test <- setAlleleMatrix(test, het.deviance.threshold = 0.1)
test <- calAlleleBoundaries(test)
test@HMM_region
