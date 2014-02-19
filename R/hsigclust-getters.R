#' @title getters for the hsigclust class
#' 
#' @description 
#' Defining getters (no setters) for the hsigclust class.
#' 
#' @aliases asimcindex, clusterList, controlFWER, hc, icovest, meigval, mpval,
#' msimeigval, nsim, plot, print, raw.data, show, summary, vsimbackvar, xmcindex
#' @genericMethods 
#' @rdname hsigclust-getters
#' @name hsigclust-getters
#' @author Patrick Kimes


#GETTERS for hsigclust class

#' @name raw.data
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return raw.data: the original data used in the analysis
setGeneric("raw.data", function(x) standardGeneric("raw.data"))
setMethod("raw.data", "hsigclust", function(x) return(x@raw.data))

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

#' @name icovest
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return icovest: the covariance estimation procedure used for the analaysis
setGeneric("icovest", function(x) standardGeneric("icovest"))
setMethod("icovest", "hsigclust", function(x) return(x@icovest))

#' @name nsim
#' @export
#' @docType methods
#' @rdname hsigclust-getters
#' @return nsim: number of simulated null datasets for each test
setGeneric("nsim", function(x) standardGeneric("nsim"))
setMethod("nsim", "hsigclust", function(x) return(x@nsim))

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



# #older approach?
# mpval <- function(x) {UseMethod('mpval', x)}
# mpval.hsigclust <- function(hsigclust) {return(hsigclust@mpval)}
