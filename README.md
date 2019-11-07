# rnaSeqNormalization
### An evolving one-stop-shop for Cory's RNA-seq preparations.

##### See vignettes (aka "articles") for working examples.

To install this package:

```
R> library(devtools)
R> biocGet("DESeq")
R> install_github("priceLab/rnaSeqNormalizer")
```

We will extend and enhance this package to accomodate a growing number
of the formats in which RNA-seq data is avaialble. 

At present (7 nov 2019) this package provides three normalization
algorithms for RNA-seq counts in a data.frame whose first column is
__ensembl_id__, and all subsequent columnns are read counts.  

Details and examples are provided in the vignettes.





