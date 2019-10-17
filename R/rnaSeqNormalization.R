#' @importFrom methods new
#'
#' @title rnaSeqNormalization
#------------------------------------------------------------------------------------------------------------------------
#' @name rnaSeqNormalization-class
#' @rdname rnaSeqNormalization-class
#' @aliases rnaSeqNormalization
#'

.rnaSeqNormalization <- setClass ("rnaSeqNormalization",
                          representation = representation(
                             mtx="matrix"
                             )
                            )

#------------------------------------------------------------------------------------------------------------------------
setGeneric('normalize', signature='obj', function(obj) standardGeneric('normalize'))
setGeneric('getNormalizedMatrix',  signature='obj', function(obj) standardGeneric('getNormalizedMatrix'))
#------------------------------------------------------------------------------------------------------------------------
#' Define an object of class rnaSeqNormalization
#'
#' @description
#' Cory, Michael and Max brewed up this method, here packaged up for easy reuse.
#'
#' @rdname rnaSeqNormalization-class
#'
#' @param mtx a matix in "raw" form
#'
#' @return A normalized matrix with dimensions, row and column names preserved
#'
#' @export
#'
#'
rnaSeqNormalization <- function(mtx)
{
   stopifnot(is.matrix(mtx))
   stopifnot(all(mtx >= 0))

   .rnaSeqNormalization(mtx=mtx)

} # ctor
#------------------------------------------------------------------------------------------------------------------------
#' get a summary of the objects
#'
#' @rdname show
#' @aliases show
#'
#' @param obj An object of class rnaSeqNormalization
#'
#' @export

setMethod('show', 'rnaSeqNormalization',

    function(object) {
       cat(sprintf ('--- rnaSeqNormalization'), '\n', sep='')
       cat(sprintf("matrix dimensions: %d x %d", nrow(object@mtx), ncol(object@mtx)))
       cat("\n")
       })

#------------------------------------------------------------------------------------------------------------------------
#' do the normalization
#'
#' @rdname normalize
#' @aliases normalize
#'
#' @param obj An object of class rnaSeqNormalization
#'
#' @export

setMethod('normalize', 'rnaSeqNormalization',

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
#' @param obj An object of class rnaSeqNormalization
#'
#' @export

setMethod('getNormalizedMatrix', 'rnaSeqNormalization',

   function(obj) {
      invisible(obj@mtx)
      })

#------------------------------------------------------------------------------------------------------------------------
