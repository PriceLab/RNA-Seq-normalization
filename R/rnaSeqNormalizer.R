#' @importFrom methods new
#' @importFrom AnnotationDbi keys select
#' @import org.Hs.eg.db
#' @import EnsDb.Hsapiens.v79
#'
#' @title rnaSeqNormalizer
#------------------------------------------------------------------------------------------------------------------------
#' @name rnaSeqNormalizer-class
#' @rdname rnaSeqNormalizer-class
#' @aliases rnaSeqNormalizer
#'

.rnaSeqNormalizer <- setClass ("rnaSeqNormalizer",
                          representation = representation(
                             mtx="matrix",
                             tbl="data.frame"
                             )
                            )

#------------------------------------------------------------------------------------------------------------------------
setGeneric('normalize', signature='obj', function(obj) standardGeneric('normalize'))
setGeneric('getNormalizedMatrix', signature='obj', function(obj) standardGeneric('getNormalizedMatrix'))
setGeneric('toGeneSymbolMatrix', signature='obj', function(obj, method) standardGeneric('toGeneSymbolMatrix'))
#------------------------------------------------------------------------------------------------------------------------
#' Define an object of class rnaSeqNormalizer
#'
#' @description
#' Cory, Michael and Max brewed up this method, here packaged up for easy reuse.
#'
#' @rdname rnaSeqNormalizer-class
#'
#' @param mtx a matix in "raw" form
#'
#' @return A normalized matrix with dimensions, row and column names preserved
#'
#' @export
#'
#'
rnaSeqNormalizer <- function(mtx)
{
   stopifnot(is.matrix(mtx))
   stopifnot(all(mtx >= 0))

   .rnaSeqNormalizer(mtx=mtx, tbl=data.frame())

} # ctor
#------------------------------------------------------------------------------------------------------------------------
#' create object of class rnaSeqNormalizer, initializing with a data.frame
#'
#'
#' @rdname data.frame.rnaSeqNormalizer
#'
#' @param mtx a matix in "raw" form
#'
#' @return A normalized matrix with dimensions, row and column names preserved
#'
#' @export
#'
#'
data.frame.rnaSeqNormalizer <- function(tbl)
{
   stopifnot(is.data.frame(tbl))

   .rnaSeqNormalizer(mtx=matrix(), tbl=tbl)

} # ctor
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
       cat(sprintf("matrix dimensions: %d x %d", nrow(object@mtx), ncol(object@mtx)))
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

setMethod('toGeneSymbolMatrix', 'rnaSeqNormalizer',

   function(obj, method) {
      browser()
      tbl.class <- .categorizeDataFrame(obj@tbl)
      switch(tbl.class,
          "ensembl column 1": return(.transformEnsemblColumn1DataFrameToGeneSymbolMatrix(obj@tbl)),
          stop(sprintf("no support yet for data frame of type '%s'"))
          )
      #mtx <- obj@mtx
      #minValue <- min(mtx[mtx > 0])
      #if(minValue == 0)
      #   minValue <- .Machine$double.eps
      #mtx.1 <- mtx + minValue
      #mtx.2 <- log10(mtx.1)
      #mtx.3 <- t(scale(t(mtx.2)))
      #means <- apply(mtx.2, 2, mean)
      #invisible(mtx.3)
      })

#------------------------------------------------------------------------------------------------------------------------
# are gene identifiers in the rownames?  in the first column?  are they gene symbols, ensembl ids?
# the goal is to
.categorizeDataFrame <- function(tbl)
{
     # create an id map data.frame to assocate ensembl gene ids with gene symbols

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
.transformEnsemblColumn1DataFrameToGeneSymbolMatrix <- function(tbl)
{
   printf("--- .transformEnsemblColumn1DataFrameToGeneSymbolMatrix")

   dim(tbl)
   mean <- apply(tbl[, -1], 1, mean)
   median <- apply(tbl[, -1], 1, median)
   sd <- apply(tbl[, -1], 1, sd)
   tbl$mean <- mean
   tbl$median <- median
   tbl$sd <- sd

   browser()
   ensg <- tbl[,1]
   tbl.map <- select(EnsDb.Hsapiens.v79, key=ensg,  columns=c("SYMBOL"),  keytype="GENEID")
   nas <- which(is.na(tbl.map$SYMBOL))
   stopifnot(length(nas) == 0)
   indices <- match(tbl[,1], tbl.map$GENEID)
   length(indices)
   head(indices)
   tbl$sym <- tbl.map$SYMBOL[indices]
   dup.syms <- which(duplicated(tbl$sym))
   length(dup.syms)
   new.order <- order(tbl$sym, tbl$sd, decreasing=TRUE)
   tbl <- tbl[new.order,]

      # discard extra rows added because a single ENSEMBL id is mapped to multiple gene symbols
   deleters <- which(duplicated(tbl$sym))
   length(deleters)
   if(length(deleters) > 0)
     tbl <- tbl[-deleters,]
   rownames(tbl) <- tbl$sym

   columns.to.delete <- c(1, match(c("mean", "median", "sd", "sym"), colnames(tbl)))
   tbl <- tbl[, -columns.to.delete]
   mtx <- as.matrix(tbl)
   return(mtx)

} # .transformEnsemblColumn1DataFrameToGeneSymbolMatrix
#------------------------------------------------------------------------------------------------------------------------
#' do the normalizer
#'
#' @rdname normalize
#' @aliases normalize
#'
#' @param obj An object of class rnaSeqNormalizer
#'
#' @export

setMethod('normalize', 'rnaSeqNormalizer',

   function(obj) {
      mtx <- obj@mtx
      minValue <- min(mtx[mtx > 0])
      if(minValue == 0)
         minValue <- .Machine$double.eps
      mtx.1 <- mtx + minValue
      mtx.2 <- log10(mtx.1)
      mtx.3 <- t(scale(t(mtx.2)))
      means <- apply(mtx.2, 2, mean)
      invisible(mtx.3)
      })

#------------------------------------------------------------------------------------------------------------------------
#' get the normalized matrix
#'
#' @rdname getNormalizedMatrix
#' @aliases getNormalizedMatrix
#'
#' @param obj An object of class rnaSeqNormalizer
#'
#' @export

setMethod('getNormalizedMatrix', 'rnaSeqNormalizer',

   function(obj) {
      invisible(obj@mtx)
      })

#------------------------------------------------------------------------------------------------------------------------
