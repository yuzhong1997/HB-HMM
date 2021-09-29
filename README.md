# HB-HMM
Extract the candidates of SNPs caused by deletion/LOH using HMM model

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

## extract the candidates of SNPs affected by deletion/LOH using HMM
```
allele_HMM <- calAlleleBoundaries(r.sub = allele_matrix$alter_less_sc, 
                                  n.sc.sub = allele_matrix$cov_sc)
```
