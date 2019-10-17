# test_rnaSeqNormalization
#------------------------------------------------------------------------------------------------------------------------
library(rnaSeqNormalization)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_basic()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))

   mtx <- get(load(system.file(package="rnaSeqNormalization", "extdata", "mtx.mayoTcx.100x300.RData")))
   normalizer <- rnaSeqNormalization(mtx)
   checkTrue(is(normalizer) == "rnaSeqNormalization")

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
test_basic <- function()
{
   message(sprintf("--- test_basic"))

   mtx <- get(load(system.file(package="rnaSeqNormalization", "extdata", "mtx.mayoTcx.100x300.RData")))
   normalizer <- rnaSeqNormalization(mtx)
   mtx.norm <- normalize(normalizer)
   checkEqualsNumeric(fivenum(mtx.norm), c(-3.59, -0.62, 0.028, 0.66, 3.23), tol=1e-2)
   checkEquals(dim(mtx), dim(mtx.norm))
   checkEquals(rownames(mtx), rownames(mtx.norm))
   checkEquals(colnames(mtx), colnames(mtx.norm))

} # test_basic
#------------------------------------------------------------------------------------------------------------------------
explore_transformation <- function()
{
   mtx.all <- get(load(system.file(package="rnaSeqNormalization", "extdata", "mtx.mayoTcx.100x300.RData")))
   mtx <- mtx.all[1:5, 1:8]
   rownames(mtx) <- paste("gene.", 1:5, sep="")
   minValue <- min(mtx[mtx > 0])
   mtx.1 <- mtx + minValue
   mtx.2 <- log10(mtx.1)
   mtx.3 <- t(scale(t(mtx.2)))
   mtx.4 <- asinh(mtx)
   quartz(); boxplot(t(mtx),   main="0) Untransformed")
   quartz(): boxplot(t(mtx.1), main="1) No zeroes")
   quartz(): boxplot(t(mtx.2), main="2) log10(mtx)")
   quartz(): boxplot(t(mtx.3), main="3) scaled & centered")
   quartz(): boxplot(t(mtx.4), main="4) asinh of untransformed")





} # explore_transformation
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
