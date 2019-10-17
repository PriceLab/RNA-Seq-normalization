# rnaSeqNormalization
### An evolving one-stop-shop for Cory's RNA-seq preparations.

##### See vignettes (aka "articles") for working examples.

To install this package:

```
R> library(devtools)
R> install_github("priceLab/rnaSeqNormalizer")
```

We will extend and enhance this package to accomodate a growing number
of the formats in which RNA-seq data is avaialble.  Identifier
transformation - for instance, from ensembl ENSG to geneSymbol - will
be supported as well.

At present (17 oct 2019) we only support one format, whose structure
can be see in this fragment:

```
> mtx[1:5, 1:8]
                X11344_TCX X11316_TCX X11431_TCX X11341_TCX X11289_TCX X11327_TCX X11334_TCX X11480_TCX
ENSG00000000003    10.6908    15.2373    13.1779    11.0598    10.4162    11.7756     7.9013     7.9966
ENSG00000000005     0.1249     0.0612     0.1011     0.1372     0.0930     0.0658     0.1624     0.1514
ENSG00000000419    23.0052    17.5626    17.0585    18.9474    16.8334    13.4531    16.9932    17.4059
ENSG00000000457    10.8656     7.2515     8.8324    10.4254     7.0682     3.5853     6.0613     7.1642
ENSG00000000460     6.8691     6.4865     5.5177     6.7045     5.6421     3.6840     4.7354     3.7082
```



