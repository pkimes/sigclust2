#' return Family-Wise Error Rate (FWER) cutoffs
#' 
#' @name fwer_cutoff-generic
#' @docType methods
#' @export
#' @keywords internal
fwer_cutoff <- function(obj, ...) {
    UseMethod("fwer_cutoff", obj)
}


#' null assumption diagnostic plots
#'
#' @name diagnostic-generic
#' @docType methods
#' @export
#' @keywords internal
diagnostic <- function(obj, ...) {
    UseMethod("diagnostic", obj)
}




