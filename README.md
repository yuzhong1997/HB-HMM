# HB-HMM
Extract the candidates of SNPs caused by deletion/LOH events using HMM model

## load sample data sets
```
load("sample_data/r.rda") # alternate allele
load("sample_data/cov.sc.rda") # total coverage
```

## load the script
```
source("HMM_code.R")
```

## set allele matrix
```
allele_matrix <- setAlleleMatrix(alter_sc = r, cov_sc = cov.sc, 
                                 het.deviance.threshold = 0.1, verbose = TRUE)
```
## plot allele-based information
```
allele_plot <- plotAlleleProfile(r.sub = allele_matrix$alter_less_sc, 
                                 n.sc.sub = allele_matrix$cov_sc, 
                                 l.sub = allele_matrix$alter_less_bulk, 
                                 n.bulk.sub = allele_matrix$cov_bulk, 
                                 snps = allele_matrix$snps)
```

## extract the candidates of continuous SNPs using HMM
```
allele_HMM <- calAlleleBoundaries(r.sub = allele_matrix$alter_less_sc, 
                                  n.sc.sub = allele_matrix$cov_sc)
```
## Update (2021/10/6):

## data preprocessing
```
alt.table <- read.table("sample_data/example.alt.matrix")
alt.table <- apply(alt.table, c(1,2), function(x) ifelse(x > 0 & x < 1, 0, x)) 
mean(alt.table > 0 & alt.table < 1)
dim(alt.table)

rownames(alt.table) <- gsub("::", ":", rownames(alt.table))

tot.table <- read.table("sample_data/example.tot.matrix") %>% as.matrix()
mean(tot.table > 0 & tot.table < 1)
dim(tot.table)

rownames(tot.table) <- gsub("::", ":", rownames(tot.table))
```

## create a object for allele-based method
```
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
test <- CreateInfercnvAlleleObject(allele_counts_matrix = alt.table,
                                   coverage_counts_matrix = tot.table,
                                   ref_annotation_file = txdb)
```

## set allele matrix
```
test <- setAlleleMatrix(test, het.deviance.threshold = 0.1)
```

## extract the candidates of continuous SNPs using HMM
```
test <- calAlleleBoundaries(test)
```
## HMM output
```
GRanges object with 24 ranges and 0 metadata columns:
       seqnames              ranges strand
          <Rle>           <IRanges>  <Rle>
   [1]    chr10     812082-84256410      *
   [2]    chr11 112045923-113240224      *
   [3]    chr11 122713154-124923820      *
   [4]    chr11     703595-72110190      *
   [5]    chr14   64442127-99909370      *
   ...      ...                 ...    ...
  [20]     chr5  10239166-109383450      *
  [21]     chr5     766098-80142489      *
  [22]     chr7 100038344-104419676      *
  [23]     chr7   1609484-159122104      *
  [24]     chrX    536656-155259147      *
  -------
  seqinfo: 24 sequences from an unspecified genome; no seqlengths
```

## Update (2021/10/16):
```
source("script/infercnv_allele.R")
file_name <- c("MGH36","MGH53", "MGH54")

res = list()
for(i in file_name){
  
  if(i == "MGH36"){
    ref_group_names = "Microglia/Macrophage" 
  }
  if(i %in% c("MGH53", "MGH54")){
    ref_group_names = c("Microglia/Macrophage","Oligodendrocytes (non-malignant)") 
  }
  
  test <- CreateInfercnvAlleleObject(allele_counts_matrix = paste0("processed_data/",i,".alt.dense.matrix"),
                                     coverage_counts_matrix = paste0("processed_data/",i,".tot.dense.matrix"),
                                     ref_annotation_file = "gencode_v19_gene_pos.txt",
                                     sample_annotation_file = paste0("processed_data/",i,".cell_types"),
                                     ref_group_names = ref_group_names,
                                     rowname_by = "::")
  test <- setAlleleMatrix(test, het.deviance.threshold = 0.1)
  test <- calAlleleBoundaries(test, distance_method = "Filter_threshold")
  plot_allele(test, name_to_plot = paste0(i,"_HMM_oct_seventeen_Filter_threshold.pdf"))
  test <- calAlleleBoundaries(test, distance_method = "Remove_NA")
  plot_allele(test, name_to_plot = paste0(i,"_HMM_oct_seventeen_Remove_NA.pdf"))
  test <- calAlleleBoundaries(test, distance_method = "HB")
  plot_allele(test, name_to_plot = paste0(i,"_HMM_oct_seventeen_HB.pdf"))
  
  res[[i]] <- test
}
```