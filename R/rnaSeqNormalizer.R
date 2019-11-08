#' @importFrom methods new
#' @importFrom AnnotationDbi keys select
#' @import org.Hs.eg.db
#' @import EnsDb.Hsapiens.v79
#' @import DESeq
#'
#' @title rnaSeqNormalizer
#------------------------------------------------------------------------------------------------------------------------
#' @name rnaSeqNormalizer-class
#' @rdname rnaSeqNormalizer-class
#' @aliases rnaSeqNormalizer
#'

.rnaSeqNormalizer <- setClass ("rnaSeqNormalizer",
                          representation = representation(
                             state="environment",
                             algorithm="character",
                             duplicate.selection.statistic="character",
                             raw.data.type="character"
                             )
                            )

#------------------------------------------------------------------------------------------------------------------------
setGeneric('normalize', signature='obj', function(obj) standardGeneric('normalize'))
setGeneric('getNormalizedMatrix', signature='obj', function(obj) standardGeneric('getNormalizedMatrix'))
setGeneric('standardizeToGeneSymbolMatrix', signature='obj', function(obj) standardGeneric('standardizeToGeneSymbolMatrix'))
#------------------------------------------------------------------------------------------------------------------------
#' Define an object of class rnaSeqNormalizer
#'
#' @description
#' Cory, Michael and Max brewed up their preferred method, to which I have added
#' asinh and vst (variance stabilizing normalization), here packaged up for easy reuse.
#' ensg identifers are mapped to geneSymbols. duplicated geneSymbols are eliminated in
#' favor of the one with the highest score according to the selected statistic:
#' mean, median, sd (standard deviation)
#'
#' @rdname rnaSeqNormalizer-class
#'
#' @param x a matix or data.frame of RNA-seq counts, rows are genes, columnns are mostly samples
#' @param algorithm a character string, one of "log+scale", "asinh", "vsn
#' @param duplicate.selection.statistc a character string, one of "median", "mean", "sd".  when
#' duplicate geneSymbols are supplied or obtained from name mapping, the highest scoring of the
#' duplicates, by the specified meansure, is kept; all others are dropped
#'
#' @return A normalized matrix with dimensions, row and column names preserved
#'
#' @export
#'
#'
rnaSeqNormalizer <- function(x, algorithm, duplicate.selection.statistic)
{
   stopifnot(algorithm %in% c("log+scale", "asinh", "vst"))
   stopifnot(duplicate.selection.statistic %in% c("mean", "median", "sd"))

   state <- new.env(parent=emptyenv())
   state$matrix <- matrix()
   state$tbl <- data.frame()

   if(is.matrix(x)){
     stop(sprintf("rnaSeqNormalizer currently supports only data.frames with ENSG ids in the first column"))
      }

   if(is.data.frame(x)){
      raw.data.type <- "data.frame"
      state$tbl <- x
      data.frame.type <- .categorizeDataFrame(x)
      if(!data.frame.type == "ensembl column 1")
         stop("rnaSeqNormalizer currently supports only data.frames with ENSG ids in the first column")
      }

   .rnaSeqNormalizer(state=state, algorithm=algorithm,
                     duplicate.selection.statistic=duplicate.selection.statistic,
                     raw.data.type=raw.data.type)

} # ctor: for data.frames with ensembl_id in first column
#------------------------------------------------------------------------------------------------------------------------
#' Define an object of class rnaSeqNormalizer from a GTEx matrix
#'
#' @description
#' Cory, Michael and Max brewed up their preferred method, to which I have added
#' asinh and vst (variance stabilizing normalization), here packaged up for easy reuse.
#' ensg identifers are mapped to geneSymbols. duplicated geneSymbols are eliminated in
#' favor of the one with the highest score according to the selected statistic:
#' mean, median, sd (standard deviation)
#'
#' @rdname rnaSeqNormalizer.gtex
#'
#' @param x a matix or data.frame of RNA-seq counts, rows are genes, columnns are mostly samples
#' @param algorithm a character string, one of "log+scale", "asinh", "vsn
#' @param duplicate.selection.statistc a character string, one of "median", "mean", "sd".  when
#' duplicate geneSymbols are supplied or obtained from name mapping, the highest scoring of the
#' duplicates, by the specified meansure, is kept; all others are dropped
#'
#' @return A normalized matrix with dimensions, row and column names preserved
#'
#' @export
#'
#'
rnaSeqNormalizer.gtex <- function(x, algorithm, duplicate.selection.statistic)
{
   stopifnot(algorithm %in% c("log+scale", "asinh", "vst"))
   stopifnot(duplicate.selection.statistic %in% c("mean", "median", "sd"))
   state <- new.env(parent=emptyenv())

   state$tbl <- data.frame()

   if(!is.matrix(x)){
     stop(sprintf("rnaSeqNormalizer.gtex currently supports only matrices with geneSymbol rownames"))
     }

   raw.data.type <- "matrix"

   if(any(duplicated(rownames(x)))){
      mean <- apply(x, 1, mean)
      median <- apply(x, 1, median)
      sd <- apply(x, 1, sd)
      x <- cbind(x, mean=mean, median=median, sd=sd)
      new.order <- order(rownames(x), x[, duplicate.selection.statistic], decreasing=TRUE)
      x <- x[new.order,]
        # discard extra rows with duplicated geneSymbols, keeping the biggest by selection.statistic
      dup.syms <- which(duplicated(rownames(x)))
      if(length(dup.syms) > 0){
         x <- x[-dup.syms,]
         }
      start.of.extra.columns <- grep("mean", colnames(x))
      end.of.extra.columns <- grep("sd", colnames(x))
      x <- x[, -(start.of.extra.columns:end.of.extra.columns)]
      dim(x)
      }

   state$matrix <- x

   .rnaSeqNormalizer(state=state, algorithm=algorithm,
                     duplicate.selection.statistic=duplicate.selection.statistic,
                     raw.data.type=raw.data.type)

} # ctor: for data.frames with ensembl_id in first column
#------------------------------------------------------------------------------------------------------------------------
#' get a summary of the objects
#'
#' @rdname show
#' @aliases show
#'
#' @param obj An object of class rnaSeqNormalizer
#'
#' @export

setMethod('show', 'rnaSeqNormalizer',

    function(object) {
       cat(sprintf ('--- rnaSeqNormalizer'), '\n', sep='')
       cat(sprintf("data.frame dimensions: %d x %d", nrow(object@state$tbl), ncol(object@state$tbl)))
       cat("\n")
       })

#------------------------------------------------------------------------------------------------------------------------
#' convert data.frame to a geneSymbol matrix
#'
#' @rdname toGeneSymbolMatrix
#' @aliases toGeneSymbolMatrix
#'
#' @param obj An object of class rnaSeqNormalizer
#' @param method A character string specifying statistic to use to resolve the case in which
#'        multiple ENSEMBL identifiers map to the same gene symbol: one of "median", "mean", "sd"
#'
#' @export

setMethod('standardizeToGeneSymbolMatrix', 'rnaSeqNormalizer',

   function(obj) {

      if(obj@raw.data.type == "data.frame"){
         tbl.class <- .categorizeDataFrame(obj@state$tbl)
         stopifnot(tbl.class %in% c("ensembl column 1"))
         if(tbl.class == "ensembl column 1")
            obj@state$matrix <- .transformEnsemblColumn1DataFrameToGeneSymbolMatrix(obj@state$tbl,
                                                                                    obj@duplicate.selection.statistic)
         }
      if(obj@raw.data.type == "matrix"){
         stop("no support yet for input matrices in standardizeGeneSymbolToMatrix")
         }
      })

#------------------------------------------------------------------------------------------------------------------------
# are gene identifiers in the rownames?  in the first column?  are they gene symbols, ensembl ids?
# the goal is to
.categorizeDataFrame <- function(tbl)
{
     # create an id map data.frame to assocate ensembl gene ids with gene symbols
   stopifnot("ensembl_id" %in% colnames(tbl))
   ensg <- tbl$ensembl_id

   tbl.map <- select(EnsDb.Hsapiens.v79, key=ensg,  columns=c("SYMBOL"),  keytype="GENEID")
     # are the rownames GENEID ids or geneSymbols
   potential.matches <- head(rownames(tbl))
   ensembl.rowname.matches <- length(which(!is.na(unlist(lapply(potential.matches,
                                                                function(id) match(id, tbl.map$GENEID))))))
   geneSymbol.rowname.matches <- length(which(!is.na(unlist(lapply(potential.matches,
                                                                function(id) match(id, tbl.map$SYMBOL))))))

     # or column 1?
   potential.matches <- head(tbl[,1])
   ensembl.column.1.matches <- length(which(!is.na(unlist(lapply(potential.matches,
                                                                function(id) match(id, tbl.map$GENEID))))))
   geneSymbol.column.1.matches <- length(which(!is.na(unlist(lapply(potential.matches,
                                                                function(id) match(id, tbl.map$SYMBOL))))))

     # just one of the above can be true: either ensg or symbols in rownames of column 1
   stopifnot(sum(c(ensembl.rowname.matches, geneSymbol.rowname.matches,
                   ensembl.column.1.matches, geneSymbol.column.1.matches)) == length(potential.matches))
   if(ensembl.rowname.matches > 0)
      return("ensembl rownames")
   if(geneSymbol.rowname.matches > 0)
      return("genesymbol rownames")
   if(ensembl.column.1.matches > 0)
      return("ensembl column 1")
   if(geneSymbol.column.1.matches > 0)
      return("genesymbol column 1")

   stop("failed to categorize data.frame")

   #printf("%d %d %d %d", ensembl.rowname.matches, geneSymbol.rowname.matches,
   #                      ensembl.column.1.matches, geneSymbol.column.1.matches)

} # .categorizeDataFrame
#------------------------------------------------------------------------------------------------------------------------
.transformEnsemblColumn1DataFrameToGeneSymbolMatrix <- function(tbl, duplicate.selection.statistic)
{
   #printf("--- .transformEnsemblColumn1DataFrameToGeneSymbolMatrix")

   dim(tbl)
   mean <- apply(tbl[, -1], 1, mean)
   median <- apply(tbl[, -1], 1, median)
   sd <- apply(tbl[, -1], 1, sd)
   tbl$mean <- mean
   tbl$median <- median
   tbl$sd <- sd

   ensg <- tbl[,1]
   tbl.map <- select(EnsDb.Hsapiens.v79, key=ensg,  columns=c("SYMBOL"),  keytype="GENEID")
   nas <- which(is.na(tbl.map$SYMBOL))
   stopifnot(length(nas) == 0)
   indices <- match(tbl[,1], tbl.map$GENEID)
   length(indices)
   head(indices)
   which(is.na(indices))
   tbl$sym <- tbl.map$SYMBOL[indices]
   mapping.failures <- which(is.na(tbl$sym))

   if(length(mapping.failures) > 0){
      tbl$sym[mapping.failures] <- tbl[mapping.failures,1]
      }

   dim(tbl)
   dup.syms <- which(duplicated(tbl$sym))
   length(dup.syms)
   new.order <- order(tbl$sym, tbl[, duplicate.selection.statistic], decreasing=TRUE)
   tbl <- tbl[new.order,]

      # discard extra rows added because a single ENSEMBL id is mapped to multiple gene symbols
   deleters <- which(duplicated(tbl$sym))
   if(length(deleters) > 0){
      tbl <- tbl[-deleters,]
      }
   dim(tbl)
   rownames(tbl) <- tbl$sym

   columns.to.delete <- c(1, match(c("mean", "median", "sd", "sym"), colnames(tbl)))
   tbl <- tbl[, -columns.to.delete]
   mtx <- as.matrix(tbl)
   return(mtx)

} # .transformEnsemblColumn1DataFrameToGeneSymbolMatrix
#------------------------------------------------------------------------------------------------------------------------
#' process according to algorithm and duplicate.selection.statistic
#'
#' @rdname getNormalizedMatrix
#' @aliases getNormalizedMatrix
#'
#' @param obj An object of class rnaSeqNormalizer
#'
#' @export

setMethod('getNormalizedMatrix', 'rnaSeqNormalizer',

   function(obj) {

      if(obj@raw.data.type == "data.frame")
         standardizeToGeneSymbolMatrix(obj)
      mtx <- obj@state$matrix

      if(obj@algorithm == "asinh")
         return(.asinh.normalize(mtx))
      if(obj@algorithm == "log+scale")
         return(.logAndScale.normalize(mtx))
      if(obj@algorithm == "vst")
         return(.vst.normalize(mtx))
      })


#------------------------------------------------------------------------------------------------------------------------
.logAndScale.normalize <- function(mtx)
{
   minValue <- min(mtx[mtx > 0])
   if(minValue == 0)
      minValue <- .Machine$double.eps

   mtx.1 <- mtx + minValue
   mtx.2 <- log10(mtx.1)
   mtx.3 <- t(scale(t(mtx.2)))
   means <- apply(mtx.2, 2, mean)

   mtx.3

} # .logAndScale.normalize
#------------------------------------------------------------------------------------------------------------------------
.asinh.normalize <- function(mtx)
{
   asinh(mtx)

} # .asinh.normalize
#------------------------------------------------------------------------------------------------------------------------
# from jim macdonald via alison paquqette: the variance statilizing transformation
.vst.normalize <- function(mtx)
{
   if(is.numeric(mtx[1,1])){
      warning("the vst algorithm needs a matrix of integer counts; converting...")
      mtx.i <- matrix(data=as.integer(mtx), nrow=nrow(mtx))
      colnames(mtx.i) <- colnames(mtx)
      rownames(mtx.i) <- rownames(mtx)
      mtx <- mtx.i
      }

   condition <- factor(rep("AD", ncol(mtx)))
   countdata <- newCountDataSet(mtx, condition) # DESseq
   countdata <- estimateSizeFactors(countdata)
   mtx.vst <- DESeq::getVarianceStabilizedData(estimateDispersions(countdata, method="blind"))
   mtx.vst

} # .vst.normalize
#------------------------------------------------------------------------------------------------------------------------
