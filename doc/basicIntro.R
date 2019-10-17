## ----setup, include = FALSE-------------------------------------------------------------------------------------------
options(width=120)
knitr::opts_chunk$set(
   collapse = TRUE,
   eval=interactive(),
   echo=TRUE,
   comment = "#>"
)

## ----loadLibraries,  results='hide'-----------------------------------------------------------------------------------
#  library(rnaSeqNormalizer)
#  library(RUnit)

## ----create, results='hide'-------------------------------------------------------------------------------------------
#  mtx <- get(load(system.file(package="rnaSeqNormalization", "extdata", "mtx.mayoTcx.100x300.RData")))
#  normalizer <- rnaSeqNormalization(mtx)

## ----run, results='show'----------------------------------------------------------------------------------------------
#  mtx.norm <- normalize(normalizer)
#  checkEqualsNumeric(fivenum(mtx.norm), c(-3.59, -0.62, 0.028, 0.66, 3.23), tol=1e-2)
#  checkEquals(dim(mtx), dim(mtx.norm))
#  checkEquals(rownames(mtx), rownames(mtx.norm))
#  checkEquals(colnames(mtx), colnames(mtx.norm))

## ----sessionInfo------------------------------------------------------------------------------------------------------
#  sessionInfo()

