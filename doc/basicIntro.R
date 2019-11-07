## ----setup, include = FALSE-------------------------------------------------------------------------------------------
options(width=120)
knitr::opts_chunk$set(
   collapse = TRUE,
   eval=TRUE,
   echo=TRUE,
   comment = "#>"
)

## ----loadLibraries,  echo=TRUE, results=FALSE, message=FALSE----------------------------------------------------------
library(rnaSeqNormalizer)
library(RUnit)

## ----readFile, results='show'-----------------------------------------------------------------------------------------
file <- system.file(package="rnaSeqNormalizer", "extdata", "tbl.ensg.column.16x10.tsv")
tbl.small <- read.table(file, sep="\t", as.is=TRUE)
checkEquals(dim(tbl.small), c(16, 10))
head(tbl.small)

## ----normalize, results='show'----------------------------------------------------------------------------------------
normalizer <- rnaSeqNormalizer(tbl.small, algorithm="asinh", duplicate.selection.statistic="median")
mtx.asinh.median <- getNormalizedMatrix(normalizer)
fivenum(mtx.asinh.median)
head(mtx.asinh.median)
hist(as.numeric(mtx.asinh.median))

## ----sessionInfo------------------------------------------------------------------------------------------------------
sessionInfo()

