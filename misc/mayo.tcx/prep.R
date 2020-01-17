# test_rnaSeqNormalizer
#------------------------------------------------------------------------------------------------------------------------
library(rnaSeqNormalizer)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
tbl <- read.table("corrected_data.tsv", sep="\t", as.is=TRUE, nrow=-1)
tbl$ensembl_id <- rownames(tbl)
rownames(tbl) <- NULL
normalizer <- rnaSeqNormalizer(tbl, "log+scale", "mean")
mtx <- getNormalizedMatrix(normalizer)
save(mtx, file="mayo.tcx.16969x262.covariateCorrection.log+scale.RData")
