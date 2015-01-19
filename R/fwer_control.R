#' return Family-Wise Error Rate (FWER) cutoffs
#' 
#' @name fwer_cutoff-generic
#' @docType methods
#' @export
#' @keywords internal
fwer_cutoff <- function(obj, ...) {
    UseMethod("fwer_cutoff", obj)
}



#' get FWER cutoffs for shc oject
#'
#' @param obj \code{shc} object
#' @param alpha numeric value specifying level
#' @param ... other parameters to be used by the function
#'
#' @name fwer_cutoff-shc
#' @export
#' @method fwer_cutoff shc
#' @author Patrick Kimes
fwer_cutoff.shc <- function(obj, alpha, ...) {
    fwer_cutoff(obj$idx_hc, alpha)
}



#' get FWER from idx_hc attribute of shc object 
#'
#' @param obj \code{shc} object
#' @param alpha numeric value specifying level
#' @param ... other parameters to be used by the function
#' 
#' @name fwer_cutoff-matrix
#' @method fwer_cutoff matrix
#' @author Patrick Kimes
#' @keywords internal
fwer_cutoff.matrix <- function(obj, alpha, ...) {
    alpha/(nrow(obj)+1) *
        apply(obj, 1, function(x) { length(unlist(x)) })
}
