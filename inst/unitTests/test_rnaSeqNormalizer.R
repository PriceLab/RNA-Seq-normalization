# test_rnaSeqNormalizer
#------------------------------------------------------------------------------------------------------------------------
library(rnaSeqNormalizer)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_basic.matrix()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))

   mtx <- get(load(system.file(package="rnaSeqNormalizer", "extdata", "mtx.mayoTcx.100x300.RData")))
   normalizer <- rnaSeqNormalizer(mtx)
   checkTrue(is(normalizer) == "rnaSeqNormalizer")

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
test_basic.matrix <- function()
{
   message(sprintf("--- test_basic.matrix"))

   mtx <- get(load(system.file(package="rnaSeqNormalizer", "extdata", "mtx.mayoTcx.100x300.RData")))
   normalizer <- rnaSeqNormalizer(mtx)
   mtx.norm <- normalize(normalizer)
   checkEqualsNumeric(fivenum(mtx.norm), c(-3.59, -0.62, 0.028, 0.66, 3.23), tol=1e-2)
   checkEquals(dim(mtx), dim(mtx.norm))
   checkEquals(rownames(mtx), rownames(mtx.norm))
   checkEquals(colnames(mtx), colnames(mtx.norm))

} # test_basic
#------------------------------------------------------------------------------------------------------------------------
test_basic.data.frame <- function()
{
   message(sprintf("--- test_basic"))

   file <- "~/github/rnaSeqNormalizer/inst/extdata/ensg.column.data.frame.16x10.RData"
   file <- system.file(package="rnaSeqNormalizer", "extdata", "ensg.column.data.frame.16x10.RData")
   checkTrue(file.exists(file))
   tbl.small <- get(load(file))
   checkEquals(dim(tbl.small), c(16, 10))
   normalizer <- data.frame.rnaSeqNormalizer(tbl.small)
   category <- rnaSeqNormalizer:::.categorizeDataFrame(tbl.small)
   checkEquals(category, "ensembl column 1")
   mtx <- rnaSeqNormalizer:::.transformEnsemblColumn1DataFrameToGeneSymbolMatrix(tbl.small)
   checkEquals(dim(mtx), c(12, 9))   # 4 ensgs eliminated, mapping to the same gene symbol

   #toGeneSymbolMatrix(normalizer, "mean")

   # mtx.norm <- normalize(normalizer)
   # checkEqualsNumeric(fivenum(mtx.norm), c(-3.59, -0.62, 0.028, 0.66, 3.23), tol=1e-2)
   # checkEquals(dim(mtx), dim(mtx.norm))
   # checkEquals(rownames(mtx), rownames(mtx.norm))
   # checkEquals(colnames(mtx), colnames(mtx.norm))

} # test_basic.data.frame
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
