# test_rnaSeqNormalizer
#------------------------------------------------------------------------------------------------------------------------
library(rnaSeqNormalizer)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_smallDataFrame()
   test_smallGTExMatrix()
   test_smallGTExMatrix_duplicatedRowNames()

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
   message(sprintf("--- test_smallDataFrame.asinh"))

   file <- system.file(package="rnaSeqNormalizer", "extdata", "tbl.ensg.column.16x10.tsv")
   tbl.small <- read.table(file, sep="\t", as.is=TRUE)
   checkEquals(dim(tbl.small), c(16, 10))

    #------------------------------------------------------------
    #   use asinh + median
    #------------------------------------------------------------

   normalizer <- rnaSeqNormalizer(tbl.small, algorithm="asinh", duplicate.selection.statistic="median")
   mtx.asinh.median <- getNormalizedMatrix(normalizer)
   fivenum(mtx.asinh.median)
   checkEqualsNumeric(fivenum(mtx.asinh.median), c(0.000000, 0.000000, 3.434929, 6.563849, 8.531491), tol=1e-2)

    #------------------------------------------------------------
    #   use asinh + sd
    #------------------------------------------------------------

   normalizer <- rnaSeqNormalizer(tbl.small, algorithm="asinh", duplicate.selection.statistic="sd")
   mtx.asinh.sd <- getNormalizedMatrix(normalizer)
   fivenum(mtx.asinh.sd)
   checkEqualsNumeric(fivenum(mtx.asinh.sd), c(0.000000, 0.000000, 3.434929, 6.563849, 8.531491), tol=1e-2)

    #------------------------------------------------------------
    #   use asinh + mean
    #------------------------------------------------------------

   normalizer <- rnaSeqNormalizer(tbl.small, algorithm="asinh", duplicate.selection.statistic="mean")
   mtx.asinh.mean <- getNormalizedMatrix(normalizer)
   fivenum(mtx.asinh.mean)
   checkEqualsNumeric(fivenum(mtx.asinh.mean), c(0.000000, 0.000000, 3.434929, 6.563849, 8.531491), tol=1e-2)
      # though the distribution is the same for these three, a slight different emerges with  mean
   checkTrue(abs(sum(mtx.asinh.median) - sum(mtx.asinh.mean)) > 0.4)

    #------------------------------------------------------------
    #   use log+scale + median
    #------------------------------------------------------------

   normalizer <- rnaSeqNormalizer(tbl.small, algorithm="log+scale", duplicate.selection.statistic="median")
   mtx.lps.median <- getNormalizedMatrix(normalizer)
   fivenum(mtx.lps.median)
   checkEqualsNumeric(fivenum(mtx.lps.median), c(-1.81, -0.76, -0.046,  0.75,  1.91),tolerance=1e-2)


    #------------------------------------------------------------
    #   use vst + median
    #------------------------------------------------------------
   suppressWarnings({
     normalizer <- rnaSeqNormalizer(tbl.small, algorithm="vst", duplicate.selection.statistic="median")
     mtx.vst.median <- getNormalizedMatrix(normalizer)
     checkEqualsNumeric(fivenum(mtx.vst.median), c(0.89, 0.89, 4.39, 8.42, 11.69), tolerance=1e-2)
     })

} # test_smallDataFrame.asinh
#------------------------------------------------------------------------------------------------------------------------
# tissue specific matrices from GTEx come as RNA-seq transcript counts, with geneSymbols already provided
# in matrix rownames.  this test ensures that rnaSeqNormalizer handles them properly
test_smallGTExMatrix <- function()
{
   message(sprintf("--- test_smallGTExMatrix"))
   file <- system.file(package="rnaSeqNormalizer", "extdata", "mtx.gtex.example.RData")
   checkTrue(file.exists(file))
   mtx <- get(load(file))

   normalizer <- rnaSeqNormalizer.gtex(mtx, algorithm="vst", duplicate.selection.statistic="median")
   suppressWarnings(
      mtx.vst.median <- getNormalizedMatrix(normalizer)
      )

   checkEquals(dim(mtx), c(100, 20))
   checkEquals(dim(mtx.vst.median), c(100, 20))
   checkEquals(rownames(mtx), rownames(mtx.vst.median))
   checkEquals(colnames(mtx), colnames(mtx.vst.median))

   checkEqualsNumeric(fivenum(mtx), c(0.00, 0.27, 4.95, 32.38, 597.50), , tolerance=1e-2)
   checkEqualsNumeric(fivenum(mtx.vst.median), c(3.73, 3.73, 4.58, 5.75, 8.95), tolerance=1e-2)

} # test_smallGTExMatrix
#------------------------------------------------------------------------------------------------------------------------
# tissue specific matrices from GTEx come as RNA-seq transcript counts, with geneSymbols already provided
# in matrix rownames.  this test ensures that rnaSeqNormalizer handles them properly
test_smallGTExMatrix_duplicatedRowNames <- function()
{
   message(sprintf("--- test_smallGTExMatrix_duplicatedRowNames"))

   file <- system.file(package="rnaSeqNormalizer", "extdata", "mtx.gtex.blood.dupRowNames.RData")
   checkTrue(file.exists(file))
   mtx <- get(load(file))

   normalizer <- rnaSeqNormalizer.gtex(mtx, algorithm="asinh", duplicate.selection.statistic="median")
   mtx.asinh.median <- getNormalizedMatrix(normalizer)

   checkEquals(dim(mtx), c(201, 20))
   checkEquals(dim(mtx.asinh.median), c(199, 20))
   checkEquals(sort(unique(rownames(mtx))), sort(unique(rownames(mtx.asinh.median))))
   checkEquals(colnames(mtx), colnames(mtx.asinh.median))

   checkEqualsNumeric(fivenum(mtx), c(0.00, 0.00, 0.27, 3.87, 4016.00), tolerance=1e-2)
   checkEqualsNumeric(fivenum(mtx.asinh.median), c(0.00, 0.00, 0.29, 2.07, 8.99), tolerance=1e-2)

} # test_smallGTExMatrix_duplicatedRowNames
#------------------------------------------------------------------------------------------------------------------------
demo.sage.matrix.1 <- function()
{
   printf("--- demo.sage.matrix.1")

   tbl.mt <- read.table("../extdata/sage/Mayo_TCX_all_counts_matrix.txt", sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
      # first 4 rows are summaries: tbl.mt$feature[1:6]
      #  "N_unmapped" "N_multimapping" "N_noFeature" "N_ambiguous" "ENSG00000223972.5" "ENSG00000227232.5"

   tbl.mt <- tbl.mt[-c(1:4),]

   colnames(tbl.mt)[1] <- "ensembl_id"
   tbl.mt$ensembl_id <- sub("\\.[0-9]+", "", tbl.mt$ensembl_id)

      #------------------------------------------------------------
      #  asinh and median on the first 100 rows and 50 columns
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.test, algorithm="asinh", duplicate.selection.statistic="median")
   mtx.1 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.1), c(100,49))
   checkEquals(fivenum(mtx.1), c(0.00, 0.00, 2.31, 5.25, 12.35), tolerance=1e-2)

      #------------------------------------------------------------
      #  log+scale and sd on the first 100 rows and 50 columns
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.test, algorithm="log+scale", duplicate.selection.statistic="median")
   mtx.2 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.2), c(100,49))
   checkEquals(fivenum(mtx.2), c(-4.94, -0.56, -0.14, 0.54, 6.86), tolerance=1e-2)

      #------------------------------------------------------------
      #  vst and mean on the first 100 rows and 50 columns
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.test, algorithm="vst", duplicate.selection.statistic="mean")
   mtx.3 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.3), c(100,49))
   checkEquals(fivenum(mtx.3), c(0.37, 0.37, 2.89, 6.57, 17.18), tolerance=1e-2)


      #------------------------------------------------------------
      #  asinh and median on the full table: 60725   279
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.mt, algorithm="asinh", duplicate.selection.statistic="median")
   mtx.4 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.4), c(58668, 278))
   checkEquals(fivenum(mtx.4), c(0.00, 0.00, 1.44, 5.14,16.75), tolerance=1e-2)

      #------------------------------------------------------------
      #  vst and mean on the full table: 60725   279
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.mt, algorithm="vst", duplicate.selection.statistic="mean")
   mtx.5 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.5), c(58668, 278))
   checkEquals(fivenum(mtx.5), c(0.34, 0.34, 2.08, 6.45, 24.87), tol=1e-2)

} # demo.sage.matrices
#------------------------------------------------------------------------------------------------------------------------
demo.sage.matrix.2 <- function()
{
   printf("--- demo.sage.matrix.2")

   tbl.mt <- read.table("../extdata/sage/MSSM_all_counts_matrix.txt", sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
      # first 4 rows are summaries: tbl.mt$feature[1:6]
      #  "N_unmapped" "N_multimapping" "N_noFeature" "N_ambiguous" "ENSG00000223972.5" "ENSG00000227232.5"

   dim(tbl.mt)
   tbl.mt <- tbl.mt[-c(1:4),]
   checkEquals(dim(tbl.mt), c(60725, 1027))

   colnames(tbl.mt)[1] <- "ensembl_id"
   tbl.mt$ensembl_id <- sub("\\.[0-9]+", "", tbl.mt$ensembl_id)

      #------------------------------------------------------------
      #  asinh and median on the first 100 rows and 50 columns
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.mt, algorithm="asinh", duplicate.selection.statistic="median")
   mtx.1 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.1), c(58668, 1026))
   checkEqualsNumeric(fivenum(mtx.1), c(0.00, 0.00,  1.44,  4.41, 15.03), tolerance=1e-2)

      #------------------------------------------------------------
      #  log+scale and mean on the full table
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.mt, algorithm="log+scale", duplicate.selection.statistic="mean")
   mtx.2 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.2), c(58668, 1026))
   checkEqualsNumeric(fivenum(mtx.2), c(-9.00, -0.48, -0.15, 0.55, 32.00), tolerance=1e-2)

      #------------------------------------------------------------
      #  vst and mean on the full table: 60725   279
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.mt, algorithm="vst", duplicate.selection.statistic="sd")
   mtx.3 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.3), c(58668, 1026))
   checkEqualsNumeric(fivenum(mtx.3), c(2.03,  2.03,  3.11,  5.61, 24.68), tolerance=1e2)

} # demo.sage.matrix.2
#------------------------------------------------------------------------------------------------------------------------
demo.sage.matrix.3 <- function()
{
   printf("--- demo.sage.matrix.3")

   tbl.mt <- read.table("../extdata/sage/Mayo_CBE_all_counts_matrix.txt", sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
      # first 4 rows are summaries: tbl.mt$feature[1:6]
      #  "N_unmapped" "N_multimapping" "N_noFeature" "N_ambiguous" "ENSG00000223972.5" "ENSG00000227232.5"

   dim(tbl.mt)
   tbl.mt <- tbl.mt[-c(1:4),]
   checkEquals(dim(tbl.mt), c(60725, 279))

   colnames(tbl.mt)[1] <- "ensembl_id"
   tbl.mt$ensembl_id <- sub("\\.[0-9]+", "", tbl.mt$ensembl_id)

      #------------------------------------------------------------
      #  asinh and median on the first 100 rows and 50 columns
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.mt, algorithm="asinh", duplicate.selection.statistic="median")
   mtx.1 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.1), c(58668, 278))
   checkEqualsNumeric(fivenum(mtx.1), c(0.00, 0.00,  1.44,  5.20, 16.46), tolerance=1e-2)

      #------------------------------------------------------------
      #  log+scale and mean on the full table
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.mt, algorithm="log+scale", duplicate.selection.statistic="mean")
   mtx.2 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.2), c(58668, 278))
   checkEqualsNumeric(fivenum(mtx.2), c(-14.93, -0.43, -0.10,   0.44, 16.61), tolerance=1e-2)

      #------------------------------------------------------------
      #  vst and mean on the full table: 60725   279
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.mt, algorithm="vst", duplicate.selection.statistic="sd")
   mtx.3 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.3), c(58668, 278))
   checkEqualsNumeric(fivenum(mtx.3), c(6.065, 6.065,  6.31,  7.65, 23.96), tolerance=1e-2)

} # demo.sage.matrix.3
#------------------------------------------------------------------------------------------------------------------------
demo.sage.matrix.4 <- function()
{
   printf("--- demo.sage.matrix.4")

   tbl.mt <- read.table("../extdata/sage/ROSMAP_all_counts_matrix.txt", sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
      # first 4 rows are summaries: tbl.mt$feature[1:6]
      #  "N_unmapped" "N_multimapping" "N_noFeature" "N_ambiguous" "ENSG00000223972.5" "ENSG00000227232.5"

   tbl.mt <- tbl.mt[-c(1:4),]
   checkEquals(dim(tbl.mt), c(60725, 640))

   colnames(tbl.mt)[1] <- "ensembl_id"
   tbl.mt$ensembl_id <- sub("\\.[0-9]+", "", tbl.mt$ensembl_id)

      #------------------------------------------------------------
      #  asinh and median on the first 100 rows and 50 columns
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.mt, algorithm="asinh", duplicate.selection.statistic="median")
   mtx.1 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.1), c(58668, 639))
   checkEqualsNumeric(fivenum(mtx.1), c(0.00, 0.00, 1.44, 4.68, 16.84), tolerance=1e-2)

      #------------------------------------------------------------
      #  log+scale and mean on the full table
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.mt, algorithm="log+scale", duplicate.selection.statistic="mean")
   mtx.2 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.2), c(58668, 639))
   checkEqualsNumeric(fivenum(mtx.2), c(-10.57, -0.49, -0.12,  0.48, 25.24), tolerance=1e-2)

      #------------------------------------------------------------
      #  vst and mean on the full table: 60725   279
      #------------------------------------------------------------

   x <- rnaSeqNormalizer(tbl.mt, algorithm="vst", duplicate.selection.statistic="sd")
   mtx.3 <-  getNormalizedMatrix(x)
      # should have one less column that tbl.test: the ensg column is gone, replaced by
      # best mapping of ensgs to geneSymbols, which are now the rownames
   checkEquals(dim(mtx.3), c(58668, 639))
   checkEqualsNumeric(fivenum(mtx.3), c(0.60, 0.60, 2.09, 5.82, 24.77), tolerance=1e2)

} # demo.sage.matrix.4
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


