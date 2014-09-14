
#' @describeIn shc get \code{data} slot
#' @aliases data,shc-method
setMethod("data", "shc", function(obj) return(obj@data))

#' @name meigval
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return meigval: matrix of sample eigenvalues from observed subtrees
setGeneric("meigval", function(x) standardGeneric("meigval"))
setMethod("meigval", "hsigclust", function(x) return(x@meigval))

#' @name msimeigval
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return msimeigval: matrix of eigenvalues used for simulating null subtrees
setGeneric("msimeigval", function(x) standardGeneric("msimeigval"))
setMethod("msimeigval", "hsigclust", function(x) return(x@msimeigval))

#' @name vsimbackvar
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return vsimbackvar: vector of background noise estimates for simulating 
#'         null subtrees
setGeneric("vsimbackvar", function(x) standardGeneric("vsimbackvar"))
setMethod("vsimbackvar", "hsigclust", function(x) return(x@vsimbackvar))

#' @name asimcindex
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return asimcindex: 
setGeneric("asimcindex", function(x) standardGeneric("asimcindex"))
setMethod("asimcindex", "hsigclust", function(x) return(x@asimcindex))

#' @name mpval
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return mpval: matrix of empirical p-values with ncols = number of cluster 
#'         indicies tested (currently 1)
setGeneric("mpval", function(x) standardGeneric("mpval"))
setMethod("mpval", "hsigclust", function(x) return(x@mpval))

#' @name mpvalnorm
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return mpvalnorm: matrix of empirical p-values with ncols = number of 
#'         cluster indicies tested (currently 1)
setGeneric("mpvalnorm", function(x) standardGeneric("mpvalnorm"))
setMethod("mpvalnorm", "hsigclust", function(x) return(x@mpvalnorm))

#' @name xmcindex
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return xmcindex:
setGeneric("xmcindex", function(x) standardGeneric("xmcindex"))
setMethod("xmcindex", "hsigclust", function(x) return(x@xmcindex))

#' @name clusterList
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return clusterList: 
setGeneric("clusterList", function(x) standardGeneric("clusterList"))
setMethod("clusterList", "hsigclust", function(x) return(x@clusterList))

#' @name hc
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return hc: 
setGeneric("hc", function(x) standardGeneric("hc"))
setMethod("hc", "hsigclust", function(x) return(x@hc))

#' @name FWERcutoffs
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return FWERcutoffs: 
setGeneric("FWERcutoffs", function(hsc, alpha) standardGeneric("FWERcutoffs"))
setMethod("FWERcutoffs", signature(hsc="hsigclust", alpha="numeric"), 
          function(hsc, alpha) {
            alpha * 
              apply(hsc@clusterList, 1, 
                    function(x) { length(unlist(x)) }) / 
              (nrow(hsc@clusterList)+1)
          })



# #older approach?
# mpval <- function(x) {UseMethod('mpval', x)}
# mpval.hsigclust <- function(hsigclust) {return(hsigclust@mpval)}
