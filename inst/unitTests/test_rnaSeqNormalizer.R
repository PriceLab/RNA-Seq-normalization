# test_rnaSeqNormalizer
#------------------------------------------------------------------------------------------------------------------------
library(rnaSeqNormalizer)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_smallDataFrame()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))

   file <- system.file(package="rnaSeqNormalizer", "extdata", "tbl.ensg.column.16x10.tsv")
   tbl.small <- read.table(file, sep="\t", as.is=TRUE)
   normalizer <- rnaSeqNormalizer(tbl.small, algorithm="asinh", duplicate.selection.statistic="median")
   checkTrue(is(normalizer) == "rnaSeqNormalizer")

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
# a common input format for expression counts is a tab-delimited file, with ensembl_id in the first column, and
# integer counts in all the other columns, which are sample names
# the desired output format is a matrx with
#   - unique geneSymbols as row names
#   - only sample ids for column names
#   - normalized numeric data throughout
#   - for any rows which, in the ensembl-to-geneSymbol mapping, have the same geneSymbols
#        * each is scored according to the chosen duplicate.selection.statistic parameter
#        * only the highest-scoring is kept, other discarded
#
test_smallDataFrame <- function()
{
   message(sprintf("--- test_smallDataFrame"))

   file <- system.file(package="rnaSeqNormalizer", "extdata", "tbl.ensg.column.16x10.tsv")
   tbl.small <- read.table(file, sep="\t", as.is=TRUE)
   checkEquals(dim(tbl.small), c(16, 10))

    #------------------------------------------------------------
    #   use asinh + median
    #------------------------------------------------------------

   normalizer <- rnaSeqNormalizer(tbl.small, algorithm="asinh", duplicate.selection.statistic="median")
   mtx.norm <- getNormalizedMatrix(normalizer)
   fivenum(mtx.norm)
   checkEqualsNumeric(fivenum(mtx.norm), c(0.000000, 0.000000, 3.434929, 6.563849, 8.531491), tol=1e-2)

    #------------------------------------------------------------
    #   use asinh + sd
    #------------------------------------------------------------

   normalizer <- rnaSeqNormalizer(tbl.small, algorithm="asinh", duplicate.selection.statistic="sd")
   mtx.norm <- getNormalizedMatrix(normalizer)
   fivenum(mtx.norm)
   checkEqualsNumeric(fivenum(mtx.norm), c(0.000000, 0.000000, 3.434929, 6.563849, 8.531491), tol=1e-2)

} # test_smallDataFrame
#------------------------------------------------------------------------------------------------------------------------
# test_basic.matrix <- function()
# {
#    message(sprintf("--- test_basic.matrix"))
#
#    mtx <- get(load(system.file(package="rnaSeqNormalizer", "extdata", "mtx.mayoTcx.100x300.RData")))
#
#
#    checkTrue(is.matrix(mtx))
#    normalizer <- rnaSeqNormalizer(mtx, algorithm="asinh", duplicate.selection.statistic="median")
#    mtx.norm <- normalize(normalizer)
#    checkEqualsNumeric(fivenum(mtx.norm), c(0.000000, 2.318597, 3.656305, 4.485816, 6.811666), tol=1e-2)
#
#    checkEquals(dim(mtx), dim(mtx.norm))
#    checkEquals(rownames(mtx), rownames(mtx.norm))
#    checkEquals(colnames(mtx), colnames(mtx.norm))
#
# } # test_basic
#------------------------------------------------------------------------------------------------------------------------
explore_transformation <- function()
{
   mtx.all <- get(load(system.file(package="rnaSeqNormalizer", "extdata", "mtx.mayoTcx.100x300.RData")))
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
create_testDataFrames <- function()
{
   tbl.big <- read.table("~/github/rnaSeqNormalizer/inst/extdata/MayoRNAseq_RNAseq_TCX_geneCounts.tsv",
                         sep="\t", header=TRUE, as.is=TRUE, nrow=-1)
   # identify the ensembl ids which
   #  - have no gene symbol
   #  - which have multiple gene symbols
   #  - are one of pairs (or triples or ...) all mapping to the same geneSymbol

   ensg <- tbl.big[,1]
   library(EnsDb.Hsapiens.v79)
   #library(EnsDb.Hsapiens.v86)

   tbl.map <- select(EnsDb.Hsapiens.v79, key=ensg,  columns=c("SYMBOL"),  keytype="GENEID")

   checkEquals(length(which(is.na(tbl.map$SYMBOL))), 0)
   checkEquals(max(table(tbl.map$GENEID)), 1)
   # checkEquals(max(table(tbl.map$SYMBOL)), 1)
   tbl.counts <- as.data.frame(table(tbl.map$SYMBOL))
   tbl.counts <- tbl.counts[order(tbl.counts$Freq, decreasing=TRUE),]
   dim(subset(tbl.counts, Freq > 1))  # 1387
   dim(subset(tbl.counts, Freq == 2)) #  818

   # look for genes with multiple ensembl ids, then identify some of those
   # in which the gene expression differs: 1162, 268
   # geneSyms.mult[c(268, 1162)]   # [1] "RPL41" "TRPV1"
   # select(EnsDb.Hsapiens.v79, key=c("RPL41", "TRPV1"),  columns=c("GENEID"),  keytype="SYMBOL")
   #   SYMBOL          GENEID
   # 1  RPL41 ENSG00000229117
   # 2  RPL41 ENSG00000273716
   # 3  RPL41 ENSG00000274700
   # 4  RPL41 ENSG00000279483
   # 5  TRPV1 ENSG00000196689
   # 6  TRPV1 ENSG00000262304

   hard.cases <- as.character(select(EnsDb.Hsapiens.v79, key=c("RPL41", "TRPV1"),  columns=c("GENEID"),  keytype="SYMBOL")$GENEID)


   #geneSyms.mult <- subset(tbl.counts, Freq > 1 & Freq < 8)$Var1
   #x <- lapply(geneSyms.mult, function(gene){printf("------ %s", gene); subset(tbl.big, ensembl_id %in% subset(tbl.map, SYMBOL==gene)$GENEID)})
      #
      # by inspection, SALL3 has three ensg ids:
   #tbl.sall3 <- select(EnsDb.Hsapiens.v79, key="SALL3",  columns=c("GENEID"),  keytype="SYMBOL")
      # all three have expression:
      #  intersect(tbl.big[,1], tbl.sall3$GENEID) # [1] "ENSG00000256463" "ENSG00000263310" "ENSG00000277015"

   #vec1 <- as.numeric(subset(tbl.big, ensembl_id=="ENSG00000256463")[, -1])   # sum: 149109
   #vec2 <- as.numeric(subset(tbl.big, ensembl_id=="ENSG00000263310")[, -1])   # sum: 0
   #vec3 <- as.numeric(subset(tbl.big, ensembl_id=="ENSG00000277015")[, -1])   # sum: 0
   #boxplot(vec1, vec2, vec3)

   set.seed(17)
   indices <- sample(seq_len(nrow(tbl.big)), 10)
   indices <- sort(indices)

   tbl.small <- tbl.big[indices,]
   dim(tbl.small)  # 10 277

   tbl.hardCases <- subset(tbl.big, ensembl_id %in% hard.cases)
   dim(tbl.hardCases)

   tbl.small <- rbind(tbl.small, tbl.hardCases)[, 1:10]
   dim(tbl.small)

   save(tbl.small, file="../extdata/ensg.column.data.frame.16x277.RData")


} # create_testDataFrames
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()


